"""Data-flow graphs are graphs where the nodes are trivial computations like
files or non-trivial computations like shell commands, and the edges are data
flows from child to parent computations. Formally, data-flow graphs are
labeled, directed, and acyclic graphs of the type
``DataFlowGraph = LabeledDirectedAcyclicGraph[Computation, ComputationLabel, DataFlowLabel]``,
where a computation label is essentially an ``n``-tuple of Python identifiers
as strings and a data-flow label is essentially a Python identifier as string.
A data flow is an edge from a computation to another computation, where the
data flows from child to parent computation. Formally, it is of the type
``DataFlow = Tuple[Computation, Computation]``.
"""

import os
import re
import shutil
import subprocess

from abc import ABCMeta, abstractmethod
from collections import Counter
from enum import Enum, unique
from inspect import signature
from multiprocessing import Process
from mypy_extensions import NamedArg
from pathlib import Path
from subprocess import PIPE, DEVNULL, CalledProcessError
from tempfile import mkdtemp
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    KeysView,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    cast,
)

from labeled_directed_acyclic_graph import LabeledDirectedAcyclicGraph, Equatable


@unique
class TransmissionKind(Enum):
    """Enumeration of possible ways to pass data between computations that run
    in separate processes.

    The less efficient but more versatile is writing and reading data to or
    from a file. It's less efficient because files are stored on the hard drive
    which has much higher access times than random access memory. It's more
    versatile because you can read files as often as you wish and even jump
    back and forth, also known as seeking, within a file when reading it. In
    other words, the data of a file is permanent.

    The more efficient but less versatile is writing and reading data to or
    from a named pipe. It's more efficient because named pipes are stored in
    random access memory. It's less versatile because you can read the data
    flowing through a pipe only once and you cannot jump back and forth when
    reading it. In other words, the data flowing through a named pipe is
    ephemeral.

    Files can be used wherever input as named pipe is possible but not the
    other way around.

    For further details see
    https://unix.stackexchange.com/questions/337739/what-is-the-difference-between-cat-file-binary-and-binary-file/337750#337750

    If you work with (named) pipes it's also helpful to understand
    https://askubuntu.com/questions/1120556/how-input-redirection-works/1120566#1120566
    and
    https://askubuntu.com/questions/172982/what-is-the-difference-between-redirection-and-pipe/1074550#1074550
    Side note: As shown in the second to last link, you can debug system calls
    using ``strace``.

    Examples:
        >>> TransmissionKind.FILE
        <TransmissionKind.FILE: 'file'>

        >>> TransmissionKind.NAMED_PIPE
        <TransmissionKind.NAMED_PIPE: 'named_pipe'>
    """

    FILE = "file"
    NAMED_PIPE = "named_pipe"


@unique
class DataFormat(Enum):
    """Enumeration of possible data formats, where values are the respective
    file extensions.

    For the list of radiance formats see
    https://floyd.lbl.gov/radiance/refer/filefmts.pdf

    Examples:
        >>> DataFormat.HDR
        <DataFormat.HDR: 'hdr'>
    """

    AMB = "amb"
    CAL = "cal"
    DAT = "dat"
    DMX = "dmx"
    FNT = "fnt"
    HDR = "hdr"
    ILL = "ill"
    MAT = "mat"
    MTX = "mtx"
    OCT = "oct"
    PIC = "pic"
    PTS = "pts"
    RAD = "rad"
    SKY = "sky"
    SMX = "smx"
    TXT = "txt"
    VEC = "vec"
    VF = "vf"  # pylint: disable=invalid-name
    WEA = "wea"
    XML = "xml"
    ZBF = "zbf"

    @classmethod
    def has_extension(cls, extension: str) -> bool:
        """Checks whether the given extension is one of a data format."""
        for data_format in DataFormat:
            if data_format.value == extension:
                return True
        return False

    @property
    def encoding(self) -> str:
        """The data format's encoding.

        If it is 'binary', you can read the corresponding file using
        ``pathlib.Path#read_bytes`` or
        ``pathlib.Path#read_text(encoding="latin-1")``. The latter works
        because the ``ascii`` compatible ``latin-1`` encoding has a (latin)
        character for each byte. For further details see
        http://python-notes.curiousefficiency.org/en/latest/python3/txt_file_processing.html#unicode-basics
        and
        http://python-notes.curiousefficiency.org/en/latest/python3/txt_file_processing.html#files-in-an-ascii-compatible-encoding-best-effort-is-acceptable

        If the encoding is not 'binary' but say ``x``, it is text and you can
        read the corresponding file using
        ``pathlib.Path#read_text(encoding=x)``. All possible encodings are
        listed on
        https://docs.python.org/3/library/codecs.html#standard-encodings


        For further details on ``pathlib`` see
        https://docs.python.org/3/library/pathlib.html
        and
        on processing files see
        http://python-notes.curiousefficiency.org/en/latest/python3/txt_file_processing.html#processing-text-files-in-python-3

        Examples:
            >>> DataFormat.HDR.encoding
            'binary'

            >>> DataFormat.CAL.encoding
            'ascii'
        """
        # pylint: disable=too-many-return-statements
        if self == DataFormat.AMB:
            return "binary"
        if self == DataFormat.CAL:
            return "ascii"
        if self == DataFormat.DAT:
            return "ascii"
        if self == DataFormat.DMX:
            return "binary"
        if self == DataFormat.FNT:
            return "ascii"
        if self == DataFormat.HDR:
            return "binary"
        if self == DataFormat.ILL:
            return "ascii"
        if self == DataFormat.MAT:
            return "ascii"
        if self == DataFormat.MTX:
            return "binary"
        if self == DataFormat.OCT:
            return "binary"
        if self == DataFormat.PIC:
            return "binary"
        if self == DataFormat.PTS:
            return "ascii"
        if self == DataFormat.RAD:
            return "ascii"
        if self == DataFormat.SKY:
            return "ascii"
        if self == DataFormat.SMX:
            return "ascii"
        if self == DataFormat.TXT:
            return "ascii"
        if self == DataFormat.VEC:
            return "ascii"
        if self == DataFormat.VF:
            return "ascii"
        if self == DataFormat.WEA:
            return "ascii"
        if self == DataFormat.XML:
            return "ascii"
        if self == DataFormat.ZBF:
            return "binary"
        raise NotImplementedError

    @property
    def is_binary(self) -> bool:
        """Whether the encoding is 'binary'.

        Examples:
            >>> DataFormat.HDR.is_binary
            True

            >>> DataFormat.CAL.is_binary
            False
        """
        return self.encoding == "binary"

    @property
    def is_text(self) -> bool:
        """Whether the encoding is not 'binary'.

        Examples:
            >>> DataFormat.HDR.is_text
            False

            >>> DataFormat.CAL.is_text
            True
        """
        return self.encoding != "binary"


class Series:
    """Series input or output meta data."""

    def __init__(self, *, format: str, wildcard: Optional[str] = None) -> None:
        """The constructor initializes series input or output meta data.

        For example, the shell command
        ```
        vwrays -vf views/south.vf -x 400 -y 400 -pj 0.7 -c 9 -ff \\
          | rfluxmtx -v -ffc -i `vwrays -vf views/south.vf -x 400 -y 400 -d` \\
          -o matrices/vmtxd/hdrIllum/south%03d.hdr -ab 1 -ad 1000 -lw 1e-4 -c 9 -n 16 \\
          - objects/GlazingVmtx.rad -i octrees/room3phDirect.oct
        ```
        outputs the series ``matrices/vmtxd/hdrIllum/south%03d.hdr`` with format
        ``%03d``. And the bash code
        ```
        for idx in {00..145}
        do
            pcomb -h -e 'ro=ri(1)*ri(2);go=gi(1)*gi(2);bo=bi(1)*bi(2)' \\
              -o matrices/vmtxd/materialMapSouth.hdr \\
              -o matrices/vmtxd/hdrIllum/south${idx}.hdr \\
              > matrices/vmtxd/hdrLum/south${idx}.hdr
        done
        ```
        cycles through the series using the wildcard ``${idx}`` and outputs
        another series again with the format ``%03d``.

        Args:
            format: The format of the series, for example, ``%03d`` means that
                the series' files are numbered ``000``, ``001``, ``002``, and so
                forth.
            wildcard: The wildcard to use instead of the given format,
                for example, ``${idx}``. This is needed, for example, when the
                actual series format is ``%03d`` but a shell command wants the
                series path to contain ``${idx}`` instead of ``%03d`` because it
                uses a loop with index variable ``idx`` to cycle through all
                files of the series. If ``None``, the given format is used as
                wildcard.

        Examples:
            >>> output_series = Series(format='%03d')
            >>> output_series.format
            '%03d'
            >>> output_series.wildcard
            '%03d'

            >>> input_series = Series(format='%03d', wildcard='${idx}')
            >>> input_series.format
            '%03d'
            >>> input_series.wildcard
            '${idx}'
        """
        # pylint: disable=redefined-builtin
        self.format = format
        self.wildcard = format if wildcard is None else wildcard

    # -----------
    # Inspired by https://stackoverflow.com/questions/45164691/recommended-way-to-implement-eq-and-hash/45170549#45170549
    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return (self.format, self.wildcard)

    def __eq__(self, other: object) -> bool:
        if type(other) is type(self):
            # pylint: disable=protected-access
            return self._members() == cast(Series, other)._members()
        return False

    def __hash__(self) -> int:
        return hash(self._members())

    # -----------

    def __str__(self) -> str:
        return "Series(format={format}, wildcard={wildcard})".format(
            format=self.format, wildcard=self.wildcard
        )

    def can_be_used_as_input_for(self, series: "Series") -> bool:
        """Whether a series with this meta data can be used as input where a
        series with the meta data ``series`` is expected, which is the case if
        and only if the two formats are identical.

        Examples:
            >>> output_series = Series(format='%03d')
            >>> input_series = Series(format='%03d', wildcard='${idx}')
            >>> output_series.can_be_used_as_input_for(input_series)
            True

            >>> output_series = Series(format='%01d')
            >>> input_series = Series(format='%05d', wildcard='${idx}')
            >>> output_series.can_be_used_as_input_for(input_series)
            False
        """
        return self.format == series.format


class InputOrOutput(metaclass=ABCMeta):
    """Abstract input or output meta data."""

    def __init__(
        self, *, data_format: DataFormat, series: Optional[Series] = None
    ) -> None:
        """The constructor initializes input or output meta data.

        Args:
            data_format: The data format of the input or output.
            series: ``None`` if the input or output is __not__ a
                series; the series meta data otherwise.
        """
        self.data_format = data_format
        self.series = series

    # -----------
    # Inspired by https://stackoverflow.com/questions/45164691/recommended-way-to-implement-eq-and-hash/45170549#45170549
    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return (self.data_format, self.series)

    def __eq__(self, other: object) -> bool:
        if type(other) is type(self):
            # pylint: disable=protected-access
            return self._members() == cast(InputOrOutput, other)._members()
        return False

    def __hash__(self) -> int:
        return hash(self._members())

    # -----------

    def __str__(self) -> str:
        return "{name}(data_format: {data_format}, series: {series})".format(
            name=self.__class__.__name__,
            data_format=self.data_format,
            series=self.series,
        )

    def file_extension(self) -> str:
        """The input's or output's file extension"""
        return cast(str, self.data_format.value)


class Input(InputOrOutput):
    """Input meta data."""

    @staticmethod
    def from_output(
        output: "Output",
        *,
        transmission_kind: Optional[TransmissionKind] = None,
        repetition_count: int = 1
    ) -> "Input":
        """Creates input meta data that is compatible with the given output
        meta data.

        Args:
            output: The output meta data that is the basis for the returned
                input meta data.
            transmission_kind: The transmission kind of the returned input meta
                data. If ``None``, the least restrictive transmission kind that
                is compatible with ``output`` is used. By default ``None``.
            repetition_count: The repetition count of the returned input meta
                data. By default ``1``.

        Returns:
            Input meta data that is comptaible with the given output meta data.

        Raises:
            AssertionError: As stated for the constructor of ``Input``.
        """
        if transmission_kind is None:
            not_none_transmission_kind = (
                TransmissionKind.NAMED_PIPE
                if output.series is None
                else TransmissionKind.FILE
            )
        else:
            not_none_transmission_kind = transmission_kind
        return Input(
            data_format=output.data_format,
            series=output.series,
            transmission_kind=not_none_transmission_kind,
            repetition_count=repetition_count,
        )

    def __init__(
        self,
        *,
        data_format: DataFormat,
        series: Optional[Series] = None,
        transmission_kind: TransmissionKind = TransmissionKind.NAMED_PIPE,
        repetition_count: int = 1
    ) -> None:
        """The constructor initializes input meta data.

        When using named pipes for transmission, to determine the required
        repetition count, you may use the helper function
        ``Helper.determine_repetition_counts``. If the required repetition count
        depends on the data or the input is a series, then you have to use
        files for transmission. And, if the repetition count is rather big, say
        100, or more than 1 and the expected input data is large, then you
        should use files also.

        Args:
            data_format: The data format of the input.
            series: ``None`` if the input is __not__ a series; the
                series meta data otherwise.
            transmission_kind: The most efficient (and least restrictive) way
                by which input data can be passed. Note that series input
                data must be passed as files.
            repetition_count: The number of times the input is read. If the
                input is passed with a file, this is of no importance and must
                be 1. However, if the input flows through a named pipe, it must
                be repeated the number of times it is read, and thus this
                number must be known to the data-flow graph executor.

        Raises:
            AssertionError: If the input is a series and the specified kind
                of transmission is __not__ ``FILE``.
            AssertionError: If the repetition count is less than or equal to 0.
            AssertionError: If the transmission kind is ``FILE`` and the
                repetition count is not 1.

        Examples:
            >>> input = Input(data_format=DataFormat.VEC)
            >>> input.data_format
            <DataFormat.VEC: 'vec'>
            >>> input.series is None
            True
            >>> input.transmission_kind
            <TransmissionKind.NAMED_PIPE: 'named_pipe'>
            >>> input.repetition_count
            1

            >>> input = Input(
            ...     data_format=DataFormat.TXT,
            ...     repetition_count=3
            ... )
            >>> input.data_format
            <DataFormat.TXT: 'txt'>
            >>> input.series is None
            True
            >>> input.transmission_kind
            <TransmissionKind.NAMED_PIPE: 'named_pipe'>
            >>> input.repetition_count
            3

            >>> series = Series(format='%03d', wildcard='${idx}')
            >>> input = Input(
            ...     data_format=DataFormat.HDR,
            ...     series=series,
            ...     transmission_kind=TransmissionKind.FILE,
            ... )
            >>> input.data_format
            <DataFormat.HDR: 'hdr'>
            >>> input.series == series
            True
            >>> input.transmission_kind
            <TransmissionKind.FILE: 'file'>
            >>> input.repetition_count
            1
        """
        assert not (
            series is not None and transmission_kind != TransmissionKind.FILE
        ), (
            "Although this input is a series its transmission kind '"
            + str(transmission_kind)
            + "' is not '"
            + str(TransmissionKind.FILE)
            + "'"  # TODO Use a data structure where this is by design impossible.
        )
        assert repetition_count >= 1, (
            "The repetition count '" + str(repetition_count) + "' is less than 1"
        )
        assert not (
            transmission_kind == TransmissionKind.FILE and repetition_count != 1
        ), (
            "Although this input transmits data through files, its repetition count '"
            "" + str(repetition_count) + "' is not 1"
        )
        super().__init__(data_format=data_format, series=series)
        self.transmission_kind = transmission_kind
        self.repetition_count = repetition_count

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return super()._members() + (self.transmission_kind, self.repetition_count)

    def __str__(self) -> str:
        return "{name}(data_format: {data_format}, series: {series}, transmission_kind: {transmission_kind}, repetition_count: {repetition_count})".format(
            name=self.__class__.__name__,
            data_format=self.data_format,
            series=self.series,
            transmission_kind=self.transmission_kind,
            repetition_count=self.repetition_count,
        )

    def with_file_transmission(self) -> "Input":
        """Returns a copy of this input meta data where the transmission kind
        is set to ``TransmissionKind.FILE`` and the repetition count is set to
        ``1``.
        """
        return Input(
            data_format=self.data_format,
            series=self.series,
            transmission_kind=TransmissionKind.FILE,
            repetition_count=1,
        )


class Output(InputOrOutput):
    """Output meta data."""

    @staticmethod
    def from_input(input: Input) -> "Output":
        """Creates output meta data that is compatible with the given input
        meta data."""
        # pylint: disable=redefined-builtin
        return Output(data_format=input.data_format, series=input.series)

    def __init__(
        self, *, data_format: DataFormat, series: Optional[Series] = None
    ) -> None:
        """The constructor initializes output meta data.

        Args:
            data_format: The data format of the output.
            series: ``None`` if the output is __not__ a series; the
                series meta data otherwise.
        """
        # pylint: disable=useless-super-delegation
        super().__init__(data_format=data_format, series=series)

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        # pylint: disable=useless-super-delegation
        return super()._members()

    def can_be_used_as_input_for(self, input: Input) -> bool:
        """Whether an output with this meta data can be used as input where an
        input with the meta data ``input`` is expected, which is the case if
        and only if the two data formats are equal and this output's optional
        series meta data can be used as input for ``input.series``.

        Examples:
            >>> output = Output(data_format=DataFormat.MTX)
            >>> one_input = Input(data_format=DataFormat.MTX)
            >>> another_input = Input(data_format=DataFormat.HDR)
            >>> output.can_be_used_as_input_for(one_input)
            True
            >>> output.can_be_used_as_input_for(another_input)
            False
        """
        # pylint: disable=redefined-builtin
        return self.data_format == input.data_format and (
            (self.series is None and input.series is None)
            or (
                self.series is not None
                and input.series is not None
                and self.series.can_be_used_as_input_for(input.series)
            )
        )


class DataFlowLabel:
    """Data-flow labels are snake-case Python identifiers as strings other than
    the reserved words ``'output'`` and ``'wildcard'`` whose format is essentially
    English words separated by underscores; for the exact definition see
    https://docs.python.org/3/reference/lexical_analysis.html#identifiers
    """

    def __init__(self, name: str) -> None:
        """The constructor initializes a data-flow label.

        Args:
            name: A snake-case Python identifier as string other than the
                reserved words ``'output'`` and ``'wildcard'``.

        Raises:
            AssertionError: If ``name`` is not a valid Python identifier or one
                of the reserved words ``'output'`` and ``'wildcard'``.

        Examples:
            >>> label = DataFlowLabel("my_data_flow_label")
            >>> str(label)
            'my_data_flow_label'

            >>> DataFlowLabel("**invalid_python_identifier**")
            Traceback (most recent call last):
            AssertionError: The name **invalid_python_identifier** is not a valid Python identifier

            >>> DataFlowLabel("output")
            Traceback (most recent call last):
            AssertionError: The name 'output' is reserved and must not be used as data-flow label.

            >>> DataFlowLabel("wildcard")
            Traceback (most recent call last):
            AssertionError: The name 'wildcard' is reserved and must not be used as data-flow label.
        """
        assert (
            name.isidentifier()
        ), "The name {name} is not a valid Python identifier".format(name=name)
        assert (
            name != "output"
        ), "The name 'output' is reserved and must not be used as data-flow label."
        assert (
            name != "wildcard"
        ), "The name 'wildcard' is reserved and must not be used as data-flow label."
        self._name = name

    # -----------
    # Inspired by https://stackoverflow.com/questions/45164691/recommended-way-to-implement-eq-and-hash/45170549#45170549
    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return (self._name,)

    def __eq__(self, other: object) -> bool:
        if type(other) is type(self):
            # pylint: disable=protected-access
            return self._members() == cast(DataFlowLabel, other)._members()
        return False

    def __hash__(self) -> int:
        return hash(self._members())

    # -----------

    def __str__(self) -> str:
        """This label's name"""
        return self._name


# Must be a list of snake-case python identifiers like ["my_group", "my_name"]
class ComputationLabel:
    """Computation labels are n-tuples of snake-case Python identifiers as
    strings whose format is essentially English words separated by underscores;
    for the exact definition see
    https://docs.python.org/3/reference/lexical_analysis.html#identifiers
    """

    def __init__(self, components: Tuple[str, ...]) -> None:
        """The constructor initializes a computation label.

        Args:
            components: A n-tuple of snake-case Python identifier as strings.

        Examples:
            >>> label = ComputationLabel(("component_one",))
            >>> str(label)
            'component_one'

            >>> label = ComputationLabel(("component_one", "component_two",))
            >>> str(label)
            'component_one/component_two'

            >>> ComputationLabel(("**non_valid_component**",))
            Traceback (most recent call last):
            AssertionError: The component **non_valid_component** is not a valid Python identifier

            >>> ComputationLabel(("component_one", "**non_valid_component**",))
            Traceback (most recent call last):
            AssertionError: The component **non_valid_component** is not a valid Python identifier
        """
        for component in components:
            assert (
                component.isidentifier()
            ), "The component {component} is not a valid Python identifier".format(
                component=component
            )
        self._components = components

    # -----------
    # Inspired by https://stackoverflow.com/questions/45164691/recommended-way-to-implement-eq-and-hash/45170549#45170549
    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return self._components

    def __eq__(self, other: object) -> bool:
        if type(other) is type(self):
            # pylint: disable=protected-access
            return self._members() == cast(ComputationLabel, other)._members()
        return False

    def __hash__(self) -> int:
        return hash(self._members())

    # -----------

    @property
    # Called `init` in Haskell
    def group_path(self) -> Tuple[str, ...]:
        """All but the last component"""
        return self._components[:-1]

    @property
    # Called `last` in Haskell
    def name(self) -> str:
        """The last component"""
        return self._components[-1]

    def to_file_name(self) -> str:
        """Converts this label to a string that is usable as part of a file
        name."""
        return "-".join(self._components)

    def __iter__(self) -> Iterator[str]:
        """Iterates over this label's components."""
        return iter(self._components)

    def __str__(self) -> str:
        """This label's components joined by '/'"""
        return "/".join(self._components)


"""A type variable for a type that is a sub-class of ``Computation``"""  # pylint: disable=pointless-string-statement
SomeComputation = TypeVar("SomeComputation", bound="Computation")


class Computation(Equatable, metaclass=ABCMeta):
    """Abstract computation representation.

    Its concrete sub-classes are nodes of data-flow graphs and they must
    implement the abstract method ``Computation.create_callable``.
    """

    @staticmethod
    def error_reporting_run(
        command: str, *, input: Optional[bytes] = None, stdout: Optional[int] = None
    ) -> bytes:
        """Runs the given command and, if given, piping it the given input,
        and, if wanted, capturing and returning its output.

        Args:
            command: The command to run.
            input: The input to pipe to the command; may be ``None``.
            stdout: How to treat the command's output: ``None`` discards it,
                ``subprocess.DEVNULL`` redirects it to the null device
                ``/dev/null``, and ``subprocess.PIPE`` captures it, see
                https://docs.python.org/3/library/subprocess.html#subprocess.DEVNULL
                and
                https://docs.python.org/3/library/subprocess.html#subprocess.PIPE
                and
                https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate

        Returns:
            The command's output if ``stdout`` is ``subprocess.PIPE``, and ``b''``
            otherwise.

        Raises:
            subprocess.CalledProcessError: If running the command fails. The
                error contains the command ``CalledProcessError#cmd``, the
                return code ``CalledProcessError#returncode``, and the error
                message ``CalledProcessError#stderr``.

        Examples:
            Run command without capturing its output:

            >>> Computation.error_reporting_run("echo __un__captured output")
            b''

            Run command piping it input:

            >>> Computation.error_reporting_run("cat", input="hello world".encode("utf-8"))
            b''

            Run command capturing its output:

            >>> Computation.error_reporting_run("echo hello world", stdout=PIPE)
            b'hello world\\n'

            Run command piping it input and capturing its output:

            >>> Computation.error_reporting_run("cat", input="hello world".encode("utf-8"), stdout=PIPE)
            b'hello world'

            Run failing command:

            >>> Computation.error_reporting_run("failing_command")
            Traceback (most recent call last):
            subprocess.CalledProcessError: Command 'set -euo pipefail && (failing_command)' returned non-zero exit status 127.
        """
        # pylint: disable=redefined-builtin
        try:
            # The command `set -euo pipefail` makes sure that errors are not silently ignored:
            # https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
            # There are other POSIX-only ways to do so:
            # https://unix.stackexchange.com/a/207678
            completed_process = subprocess.run(
                "set -euo pipefail && ({command})".format(command=command),
                input=input,
                stdout=stdout,
                stderr=PIPE,
                check=True,
                shell=True,
                executable="/bin/bash",
            )
            if stdout is not None and stdout is not DEVNULL:
                return completed_process.stdout
            return "".encode("utf-8")
        except CalledProcessError as error:
            print(
                "The command '"
                + error.cmd
                + "' exited with return code '"
                + str(error.returncode)
                + "' and error message '"
                + str(error.stderr)
                + "'"
            )
            raise error

    @staticmethod
    def _output_paths_for_each_repetition_count(  # pylint: disable=invalid-name
        output_paths_with_repetition_counts: Set[Tuple[Path, int]]
    ) -> List[List[Path]]:
        """Collects output paths for repetition counts ``<= 1``, ``<= 2``, ... and
        so on"""
        max_repetition_count = max(
            [
                repetition_count
                for _, repetition_count in output_paths_with_repetition_counts
            ]
        )
        return [
            [
                output_path
                for output_path, repetition_count in output_paths_with_repetition_counts
                if repetition_count >= i + 1
            ]
            for i in range(max_repetition_count)
        ]

    @staticmethod
    def _repeat_output(  # pylint: disable=invalid-name
        output: bytes, output_paths_for_each_repetition_count: List[List[Path]]
    ) -> None:
        """Sends the output to the given paths the requested number of
        times."""
        for output_paths in output_paths_for_each_repetition_count:
            Computation.error_reporting_run(
                "tee " + " ".join(map(str, output_paths)), input=output, stdout=DEVNULL
            )

    def __init__(
        self, *, label_to_input: Dict[DataFlowLabel, Input], output: Optional[Output]
    ) -> None:
        """The constructor initializes a new (abstract) computation.

        Args:
            label_to_input: The map from labels to to the computation's inputs.
            output: ``None`` if the computation does __not__ have
                output; the output meta data otherwise.
        """
        self._label_to_input = label_to_input
        self._output = output

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return (self._label_to_input, self.output)

    def __str__(self) -> str:
        return "Computation(label_to_input={label_to_input}, output={output})".format(
            label_to_input={
                str(label): str(inpud) for label, inpud in self._label_to_input.items()
            },
            output=self._output,
        )

    @abstractmethod
    def with_file_inputs_for(
        self: SomeComputation, input_labels: Set[DataFlowLabel]
    ) -> SomeComputation:
        """Returns a copy of this computation where the input meta data
        corresponding to the given labels is adapted such that the transmission
        kind is set to ``TransmissionKind.FILE`` and the repetition count is
        set to ``1``.

        This method must be implemented by sub-classes. To do so, they may use
        the method ``_label_to_input_with_file_inputs_for``.
        """

    def _label_to_input_with_file_inputs_for(
        self, input_labels: Set[DataFlowLabel]
    ) -> Dict[DataFlowLabel, Input]:
        """Returns a copy of ``_label_to_input`` where the input meta data
        corresponding to the given labels is adapted such that the transmission
        kind is set to ``TransmissionKind.FILE`` and the repetition count is set
        to ``1``.
        """
        return {
            label: (inpud.with_file_transmission() if label in input_labels else inpud)
            for label, inpud in self._label_to_input.items()
        }

    @abstractmethod
    def create_callable(
        self, computation_label: ComputationLabel, *, debug: bool = False
    ) -> Callable[
        [
            NamedArg(Dict[DataFlowLabel, Path], "label_to_input_path"),
            NamedArg(Set[Tuple[Path, int]], "output_paths_with_repetition_counts"),
        ],
        None,
    ]:
        """Creates a callable for this computation that when invoked performs
        this computation. It must have two precisely two named parameters,
        namely ``label_to_input_path: Dict[DataFlowLabel, Path]`` and
        ``output_paths_with_repetition_counts: Set[Tuple[Path, int]]``.

        The first parameter is a map from input labels to input paths from
        which input can be read. The input path points either to a named pipe
        or a file depending on what the input meta data for the given label
        wanted and what the source computation was able to give. For example,
        if the input meta data says that named pipes are possible and the
        source computation outputs a stream of data, then that stream is piped
        through a named pipe whose path will be given (before or after data
        flow has begun). But, if for some reason the source outputs a file,
        then the corresponding file path will be given (after all data has been
        written to the file).

        The second parameter is a set of output paths and how often the output
        is required. The callable must write output to the output path as often
        as it is required.

        This method must be implemented by sub-classes.

        Args:
            computation_label: This computation's label within the data-flow
                graph for whose execution the returned callable is being used.
            debug: Whether to use cached results in the directory ``./tmp`` and
                print debug output.
        """

    @property
    def input_labels(self) -> KeysView[DataFlowLabel]:
        """All input labels through which this computation expects input."""
        return self._label_to_input.keys()

    def optional_input(self, label: DataFlowLabel) -> Optional[Input]:
        """Returns the input meta data for the given input label if this
        computation expects input through it and ``None`` otherwise."""
        return self._label_to_input.get(label, None)

    def has_input(self, label: DataFlowLabel) -> bool:
        """Checks whether this computation expects input through the given
        input label."""
        return label in self._label_to_input

    def input(self, label: DataFlowLabel) -> Input:
        """Returns the input meta data this computation expects through the
        given input label.

        Raises:
            AssertionError: If this computation does not expect input through
                the given input label.
        """
        assert label in self._label_to_input, (
            "There is no input with label '"
            + str(label)
            + "' in the computation '"
            + str(self)
            + "'"
        )
        return self._label_to_input[label]

    @property
    def output(self) -> Optional[Output]:
        """The output meta data if this computation has output and ``None``
        otherwise"""
        return self._output

    @property
    def meta_data(self) -> Dict[str, Any]:
        """Meta data that is neither input nor output meta data (should be
        extended by sub-classes and is only used in the visualization of the
        data-flow graph)"""
        return {}

    def is_equivalent_to(self, other: object) -> bool:
        """Checks whether this computation is equivalent to the other one.

        Sub-classes should not override this method but rather just extend the
        method ``Computation._members`` to include all members that uniquely
        identify this computation. Its result is used by this method to check
        equivalence.

        Args:
            other: The other computation.
        """
        if type(other) is type(self):
            # pylint: disable=protected-access
            return self._members() == cast(Computation, other)._members()
        return False

    # private

    def _assert_existence_of_output_paths(  # pylint: disable=invalid-name
        self, output_paths_with_repetition_counts: Set[Tuple[Path, int]]
    ) -> None:
        """Assert that the given paths exist, that is, point to a named pipe,
        file, or directory."""
        for (output_path, _) in output_paths_with_repetition_counts:
            assert output_path.exists(), (
                "There is no file or pipe with path '"
                + str(output_path)
                + "' to which the output of the computation '"
                + str(self)
                + "' can be written'"
            )


class File(Computation):
    """File representation."""

    def __init__(
        self,
        path: Path,
        *,
        data_format: Optional[DataFormat] = None,
        series: Optional[Series] = None
    ) -> None:
        """The constructor initializes a new file.

        The output, including data format and possible series format is deduced
        from the path.

        Args:
            path: The file's path.
            data_format: The data format of the file. If ``None``, the path's
                extension is used to determine the format. By default ``None``.
            series: The series meta data of the possible series of files. If
                ``None``, and if the path contains a pattern that matches the
                regular expression ``%[0-9]+d``, then it's supposed to be
                a series of files, where the format is the matching pattern;
                and otherwise it's supposed to be a good old file. If not
                ``None``, then the path must contain the named placeholder
                ``{wildcard}`` which, when the series is being used as input of
                another computation, is replaced by the wildcard that
                computation specified for the input. By default ``None``.

        Raises:
            AssertionError: If no data format is given and the file path does
                not have an extension or the extension does not belong to a
                known data format.
            AssertionError: If series meta data is given and the file path does
                not contain the named placeholder ``{wildcard}``.
            AssertionError: If the file is not a series and the path does not
                exist, or if the file is a series and its directory does not
                exist.

        Examples:
            >>> tmp = getfixture('tmp_path')

            >>> txt_path = tmp / "my.txt"
            >>> txt_path.touch()
            >>> txt_file = File(txt_path)
            >>> str(txt_file.path) # doctest: +ELLIPSIS
            '.../my.txt'
            >>> txt_file.output.data_format
            <DataFormat.TXT: 'txt'>
            >>> txt_file.output.series is None
            True

            >>> mtx_path = tmp / "my"
            >>> mtx_path.touch()
            >>> mtx_file = File(mtx_path, data_format=DataFormat.MTX)
            >>> str(mtx_file.path) # doctest: +ELLIPSIS
            '.../my'
            >>> mtx_file.output.data_format
            <DataFormat.MTX: 'mtx'>
            >>> mtx_file.output.series is None
            True

            >>> hdr_path = tmp / "hdr-series" / "my-%03d.hdr"
            >>> hdr_path.parent.mkdir()
            >>> hdr_series = File(hdr_path)
            >>> str(hdr_series.path) # doctest: +ELLIPSIS
            '.../hdr-series/my-{wildcard}.hdr'
            >>> hdr_series.output.data_format
            <DataFormat.HDR: 'hdr'>
            >>> hdr_series.output.series.format
            '%03d'
            >>> hdr_series.output.series.wildcard
            '%03d'

            >>> vec_path = tmp / "vec-series" / "my-{wildcard}.vec"
            >>> vec_path.parent.mkdir()
            >>> vec_series = File(vec_path, series=Series(format="%04d"))
            >>> str(vec_series.path) # doctest: +ELLIPSIS
            '.../vec-series/my-{wildcard}.vec'
            >>> vec_series.output.data_format
            <DataFormat.VEC: 'vec'>
            >>> vec_series.output.series.format
            '%04d'
            >>> vec_series.output.series.wildcard
            '%04d'

            >>> File(Path("no-extension"))
            Traceback (most recent call last):
            AssertionError: The file path 'no-extension' does not have an extension; it is needed to determine the file's data format.

            >>> File(Path("my.unknown-extension"))
            Traceback (most recent call last):
            AssertionError: The file path 'my.unknown-extension' has an extension, namely 'unknown-extension', that does not belong to a data format.

            >>> File(Path("no-wildcard.hdr"), series=Series(format="%04d"))
            Traceback (most recent call last):
            AssertionError: The file path 'no-wildcard.hdr' does not contain the named placeholder '{wildcard}' although the series meta data 'Series(format=%04d, wildcard=%04d)' is given.

            >>> File(Path("non-existent-file.txt"))
            Traceback (most recent call last):
            AssertionError: There is no file with path 'non-existent-file.txt'.

            >>> File(Path("non-existent/series-%03d.txt"))
            Traceback (most recent call last):
            AssertionError: There is no directory with path 'non-existent'.
        """
        if data_format is None:
            extension = path.suffix[1:]
            assert extension, (
                "The file path '"
                + str(path)
                + "' does not have an extension; it is needed to determine the file's data format."
            )
            assert DataFormat.has_extension(extension), (
                "The file path '"
                + str(path)
                + "' has an extension, namely '"
                + str(extension)
                + "', that does not belong to a data format."
            )
            data_format = DataFormat(extension)
        if series is None:
            match = re.search("%[0-9]+d", str(path))
            if match is not None:
                series = Series(format=match[0])
                path = Path(
                    str(path)[: match.start()] + "{wildcard}" + str(path)[match.end() :]
                )
        else:
            assert "{wildcard}" in str(path), (
                "The file path '"
                + str(path)
                + "' does not contain the named placeholder '{wildcard}' although the series meta data '"
                + str(series)
                + "' is given."
            )
        if series is None:
            assert path.exists(), "There is no file with path '" + str(path) + "'."
        else:
            assert path.parent.exists(), (
                "There is no directory with path '" + str(path.parent) + "'."
            )
        super().__init__(
            label_to_input={},
            output=Output(
                # The suffix `path.suffix()` includes `'.'`, for example,
                # `Path('my.txt').suffix()` is equal to `'.txt'`. We remove the
                # dot by taking the slice `[1:]`.
                data_format=data_format,
                series=cast(Series, series),
            ),
        )
        self.path = path

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return super()._members() + (self.path,)

    def __str__(self) -> str:
        return "File(path={path}) < {computation}".format(
            path=self.path, computation=super().__str__()
        )

    @property
    def meta_data(self) -> Dict[str, Any]:
        """See ``Computation.meta_data``"""
        return {**super().meta_data, **{"path": self.path}}

    def with_file_inputs_for(self, input_labels: Set[DataFlowLabel]) -> "File":
        """Returns itself"""
        return self

    def create_callable(
        self, computation_label: ComputationLabel, *, debug: bool = False
    ) -> Callable[
        [
            NamedArg(Dict[DataFlowLabel, Path], "label_to_input_path"),
            NamedArg(Set[Tuple[Path, int]], "output_paths_with_repetition_counts"),
        ],
        None,
    ]:
        """This method must not bee invoked!

        Raises:
            AssertionError: Always.
        """
        assert False


class ShellCommand(Computation):
    """Shell command representation."""

    def __init__(  # pylint: disable=invalid-name
        self,
        formatted_command: str,
        *,
        label_to_input: Dict[DataFlowLabel, Input],
        output: Optional[Output],
        redirect_inputs_to_files_and_discard_command: bool = False
    ) -> None:
        """The constructor initializes a new shell command.

        Args:
            formatted_command: The shell command to run with at least one named
                placeholder for each input, whose name is the input's label.
                And, if the command's output is a series, there can and must
                exist the named placeholder ``'output'`` for the output. The input
                and output placeholders will be replaced by named pipe or file
                paths before execution; from those the inputs are meant to be
                read and the output is meant to be written to.
            label_to_input: The map from labels to to the command's inputs.
            output: The output of the command; may be ``None``.
            redirect_inputs_to_files_and_discard_command: Whether to redirect
                inputs to files in the directory ``./tmp`` and discard the
                command. Helpful for debugging.

        Raises:
            AssertionError: If there is an input label without corresponding
                named placeholder in the formatted command.
            AssertionError: If the output of the command is a series and the
                formatted command does not contain the named placeholder
                'output'.

        Examples:
            >>> sleep = ShellCommand(
            ...     "sleep 42",
            ...     label_to_input={},
            ...     output=None,
            ... )
            >>> sleep.formatted_command
            'sleep 42'

            >>> echo = ShellCommand(
            ...     "echo hello world",
            ...     label_to_input={},
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> echo.formatted_command
            'echo hello world'

            >>> count_words = ShellCommand(
            ...     "wc --words {a}",
            ...     label_to_input={DataFlowLabel("a"): Input(data_format=DataFormat.TXT),},
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> count_words.formatted_command
            'wc --words {a}'

            >>> merge = ShellCommand(
            ...     "paste {first_list} {second_list}",
            ...     label_to_input={
            ...         DataFlowLabel("first_list"): Input(data_format=DataFormat.TXT),
            ...         DataFlowLabel("second_list"): Input(data_format=DataFormat.TXT),
            ...     },
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> merge.formatted_command
            'paste {first_list} {second_list}'
        """
        assert not (
            output is not None
            and output.series is not None
            and "{output}" not in formatted_command
        ), (
            "Although the output of this shell command is a series, the formatted command '"
            + str(formatted_command)
            + "' does not contain the named placeholder '{output}'"
        )
        for input_label in label_to_input:
            assert "{" + str(input_label) + "}" in formatted_command, (
                "The input with label '"
                + str(input_label)
                + "' is not used in the command '"
                + formatted_command
                + "'"
            )
        super().__init__(label_to_input=label_to_input, output=output)
        self.formatted_command = formatted_command
        self.redirect_inputs_to_files_and_discard_command = (
            redirect_inputs_to_files_and_discard_command
        )

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return super()._members() + (
            self.formatted_command,
            self.redirect_inputs_to_files_and_discard_command,
        )

    def __str__(self) -> str:
        return "ShellCommand(formatted_command={formatted_command}) < {computation}".format(
            formatted_command=self.formatted_command, computation=super().__str__()
        )

    @property
    def meta_data(self) -> Dict[str, Any]:
        """See ``Computation.meta_data``"""
        return {
            **super().meta_data,
            **{
                "formatted_command": self.formatted_command,
                "redirect_inputs_to_files_and_discard_command": self.redirect_inputs_to_files_and_discard_command,
            },
        }

    def with_file_inputs_for(self, input_labels: Set[DataFlowLabel]) -> "ShellCommand":
        """See ``Computation.with_file_inputs_for``"""
        if input_labels:
            return ShellCommand(
                self.formatted_command,
                label_to_input=self._label_to_input_with_file_inputs_for(input_labels),
                output=self.output,
                redirect_inputs_to_files_and_discard_command=self.redirect_inputs_to_files_and_discard_command,
            )
        return self

    def create_callable(
        self, computation_label: ComputationLabel, *, debug: bool = False
    ) -> Callable[
        [
            NamedArg(Dict[DataFlowLabel, Path], "label_to_input_path"),
            NamedArg(Set[Tuple[Path, int]], "output_paths_with_repetition_counts"),
        ],
        None,
    ]:
        """See ``Computation.create_callable``"""

        def kallable(  # pylint: disable=invalid-name
            *,
            label_to_input_path: Dict[DataFlowLabel, Path],
            output_paths_with_repetition_counts: Set[Tuple[Path, int]]
        ) -> None:
            if self.redirect_inputs_to_files_and_discard_command:
                ShellCommand._redirect_inputs_to_files(
                    computation_label, label_to_input_path
                )
                print(
                    "Done redirecting inputs to files and discarding command '"
                    + self.formatted_command
                    + "'"
                )
            elif self.output is None:
                assert not output_paths_with_repetition_counts, (
                    "Although the shell command '"
                    + str(self)
                    + "' does not have output"
                    " it has been given the output paths and repetition counts '"
                    + ", ".join(
                        [
                            str(output_path) + " times " + str(count)
                            for output_path, count in output_paths_with_repetition_counts
                        ]
                    )
                    + "'"
                )
                ShellCommand._run_command_without_output(
                    self.formatted_command, label_to_input_path=label_to_input_path
                )
            else:
                assert len(output_paths_with_repetition_counts) >= 1, (
                    "There is no output path to which the output of the formatted command '"
                    + self.formatted_command
                    + "' can be written"
                )
                # Note that using cached results results in the inputs not being consumed!
                # TODO Do not use temporary output paths. Maybe use the working directory instead?
                temporary_output_path = Path("tmp", *computation_label)
                if self.output is not None and self.output.series is not None:
                    ShellCommand._run_command_with_series_output(
                        self.formatted_command,
                        label_to_input_path=label_to_input_path,
                        output_paths_with_repetition_counts=output_paths_with_repetition_counts,
                        temporary_output_path=temporary_output_path,
                        debug=debug,
                    )
                else:
                    self._assert_existence_of_output_paths(
                        output_paths_with_repetition_counts
                    )
                    ShellCommand._run_command_with_normal_output(
                        self.formatted_command,
                        label_to_input_path=label_to_input_path,
                        output_paths_with_repetition_counts=output_paths_with_repetition_counts,
                        temporary_output_path=temporary_output_path,
                        debug=debug,
                    )

        return kallable

    @staticmethod
    def _redirect_inputs_to_files(
        computation_label: ComputationLabel,
        label_to_input_path: Dict[DataFlowLabel, Path],
    ) -> None:
        for input_label, input_path in label_to_input_path.items():
            output_path = Path("tmp", *computation_label, str(input_label))
            output_path.parent.mkdir(parents=True, exist_ok=True)
            Computation.error_reporting_run(
                "cp {inpud} {output}".format(inpud=input_path, output=output_path)
            )

    @staticmethod
    def _run_command_without_output(  # pylint: disable=invalid-name
        formatted_command: str, *, label_to_input_path: Dict[DataFlowLabel, Path]
    ) -> None:
        """Runs the given formatted command."""
        command = formatted_command.format(
            **ShellCommand._convert_label_to_input_path_map_into_str_to_path_map(
                label_to_input_path
            )
        )
        Computation.error_reporting_run(command)

    @staticmethod
    def _run_command_with_series_output(  # pylint: disable=invalid-name
        formatted_command: str,
        *,
        label_to_input_path: Dict[DataFlowLabel, Path],
        output_paths_with_repetition_counts: Set[Tuple[Path, int]],
        temporary_output_path: Path,
        debug: bool
    ) -> None:
        """Runs the given formatted command and saves the series output."""
        assert len(output_paths_with_repetition_counts) == 1, (
            "The set of output paths with repetition counts for the command '"
            + formatted_command
            + "', namely '"
            + str(output_paths_with_repetition_counts)
            + "' contains more than 1 element"
        )
        # `next(iter(set))` returns one element of a set without removing it
        output_path, _ = next(iter(output_paths_with_repetition_counts))
        if debug and temporary_output_path.exists():
            print(
                "Using cached result '"
                + str(temporary_output_path)
                + " for command "
                + formatted_command
                + "'"
            )
            Computation.error_reporting_run(
                "cp {input}/* {output}".format(
                    input=temporary_output_path, output=output_path.parent
                )
            )
        else:
            command_arguments = {
                **{"output": output_path},
                **ShellCommand._convert_label_to_input_path_map_into_str_to_path_map(
                    label_to_input_path
                ),
            }  # Merge two dictionaries
            command = formatted_command.format(
                **command_arguments  # Make dictionary key-word arguments
            )
            Computation.error_reporting_run(command)

    @staticmethod
    def _run_command_with_normal_output(  # pylint: disable=invalid-name
        formatted_command: str,
        *,
        label_to_input_path: Dict[DataFlowLabel, Path],
        output_paths_with_repetition_counts: Set[Tuple[Path, int]],
        temporary_output_path: Path,
        debug: bool
    ) -> None:
        """Runs the given formatted command and saves the file output."""
        output_paths_for_each_repetition_count = Computation._output_paths_for_each_repetition_count(
            output_paths_with_repetition_counts
        )
        if debug and temporary_output_path.is_file():
            print(
                "Using cached result '"
                + str(temporary_output_path)
                + " for command "
                + formatted_command
                + "'"
            )
            command = "cat {path}".format(path=temporary_output_path)
        else:
            command = formatted_command.format(
                **ShellCommand._convert_label_to_input_path_map_into_str_to_path_map(
                    label_to_input_path
                )
            )
        redirecting_command = "( {command} ) | tee {paths}".format(
            command=command,
            paths=" ".join(map(str, output_paths_for_each_repetition_count[0])),
        )
        stdout = PIPE if len(output_paths_for_each_repetition_count) >= 2 else DEVNULL
        output = Computation.error_reporting_run(redirecting_command, stdout=stdout)
        Computation._repeat_output(output, output_paths_for_each_repetition_count[1:])

    @staticmethod
    def _convert_label_to_input_path_map_into_str_to_path_map(
        label_to_input_path: Dict[DataFlowLabel, Path]
    ) -> Dict[str, Path]:
        return {
            str(label): input_path for label, input_path in label_to_input_path.items()
        }


class Function(Computation):
    """Python function representation."""

    @staticmethod
    def _make_input_parameter_name(label: DataFlowLabel) -> str:
        """Makes the input parameter name of the given label."""
        return "{}_input_path".format(label)

    @staticmethod
    def _assert_validity(
        function: Callable[..., Optional[bytes]],
        *,
        label_to_input: Dict[DataFlowLabel, Input],
        output: Optional[Output]
    ) -> None:
        """Asserts that the function has a corresponding named parameter for
        each input and the named parameter ``'output_path'`` if the function
        has series output."""
        expected_input_paramter_names = [
            Function._make_input_parameter_name(label)
            for label in label_to_input.keys()
        ]
        expected_output_parameter_names = (
            [] if output is None or output.series is None else ["output_path"]
        )
        expected_paramter_names = (
            expected_input_paramter_names + expected_output_parameter_names
        )
        actual_paramter_names = signature(function).parameters.keys()
        assert Counter(expected_paramter_names) == Counter(actual_paramter_names), (
            "The expected named parameters "
            + str(expected_paramter_names)
            + " are not the same as the actual named parameters "
            + str(actual_paramter_names)
            + " of the function "
            + str(function)
        )

    def __init__(
        self,
        function: Callable[..., Optional[bytes]],
        *,
        label_to_input: Dict[DataFlowLabel, Input],
        output: Optional[Output]
    ) -> None:
        """The constructor initializes a new function.

        Args:
            function: The function to invoke. For each input label, the
                function must have a corresponding parameter of type ``Path``
                whose name is equal to the label followed by ``'_input_path'``,
                for example, if there is an input label ``DataFlowLabel('my')``,
                there must be a parameter named ``'my_input_path'`` of type
                ``Path``.  If the function has non-series output, it must return
                it as a byte sequence; if you have a UTF-8 string ``x``, you can
                encode it as byte sequence with ``x.encode('utf-8')``; for
                details see
                https://docs.python.org/3/library/stdtypes.html#str.encode and
                https://docs.python.org/3/library/codecs.html#standard-encodings
                If the function has series output, it must have a named
                parameter ``'output_path'`` of type ``Path`` whose value will be
                a path containing the specified series wildcard (which defaults
                to the series format) and it must write output to files whose
                paths are compatible with the output path and series format,
                for example, if the output path is ``'./names/name-*.txt'`` and
                the format is ``'%01d'``, compatible file paths are
                ``'./names/name-0.txt'``, ``'./names/name-1.txt'``, and so on until
                ``'./names/name-9.txt'``.
            label_to_input: The map from labels to to the function's inputs.
            output: The output of the function; may be ``None``.

        Raises:
            AssertionError: If there is an input label without corresponding
                parameter.

        Examples:
            >>> print = Function(
            ...     lambda my_input_path: print(pathlib.Path(my_input_path).read_text()),
            ...     label_to_input={DataFlowLabel("my"): Input(data_format=DataFormat.VEC)},
            ...     output=None,
            ... )

            >>> def replace_function(
            ...     my_input_path: Path,
            ... ) -> Optional[bytes]:
            ...     inpud = Path(my_input_path).read_text()
            ...     return inpud.replace("my", "yours").encode("utf-8")
            >>> replace = Function(
            ...     function=replace_function,
            ...     label_to_input={
            ...         DataFlowLabel("my"): Input(data_format=DataFormat.TXT),
            ...     },
            ...     output=Output(data_format=DataFormat.MTX),
            ... )
        """
        Function._assert_validity(
            function, label_to_input=label_to_input, output=output
        )
        super().__init__(label_to_input=label_to_input, output=output)
        self.function = function

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return super()._members() + (self.function,)

    def __str__(self) -> str:
        return "Function(function={function}) < {computation}".format(
            function=self.function, computation=super().__str__()
        )

    @property
    def meta_data(self) -> Dict[str, Any]:
        """See ``Computation.meta_data``"""
        return {**super().meta_data, **{"function": self.function}}

    def with_file_inputs_for(self, input_labels: Set[DataFlowLabel]) -> "Function":
        """See ``Computation.with_file_inputs_for``"""
        if input_labels:
            return Function(
                self.function,
                label_to_input=self._label_to_input_with_file_inputs_for(input_labels),
                output=self.output,
            )
        return self

    def create_callable(
        self, computation_label: ComputationLabel, *, debug: bool = False
    ) -> Callable[
        [
            NamedArg(Dict[DataFlowLabel, Path], "label_to_input_path"),
            NamedArg(Set[Tuple[Path, int]], "output_paths_with_repetition_counts"),
        ],
        None,
    ]:
        """See ``Computation.create_callable``"""

        def kallable(  # pylint: disable=invalid-name
            *,
            label_to_input_path: Dict[DataFlowLabel, Path],
            output_paths_with_repetition_counts: Set[Tuple[Path, int]]
        ) -> None:
            input_parameter_name_to_path = {
                Function._make_input_parameter_name(label): input_path
                for label, input_path in label_to_input_path.items()
            }
            if self.output is None:
                assert not output_paths_with_repetition_counts, (
                    "Although the function '" + str(self) + "' does not have output"
                    " it has been given the output paths and repetition counts '"
                    + ", ".join(
                        [
                            str(output_path) + " times " + str(count)
                            for output_path, count in output_paths_with_repetition_counts
                        ]
                    )
                    + "'"
                )
                output = self.function(**input_parameter_name_to_path)
                assert output is None, (
                    "The function '"
                    + str(self)
                    + "' returned output although it declared not to do so"
                )
            else:
                assert output_paths_with_repetition_counts, (
                    "There is no output path to which the output of the function '"
                    + str(self)
                    + "' can be written"
                )
                if self.output.series is None:
                    self._assert_existence_of_output_paths(
                        output_paths_with_repetition_counts
                    )
                    output = self.function(**input_parameter_name_to_path)
                    assert output is not None, (
                        "The function '"
                        + str(self)
                        + "' did not return output although it declared to do so"
                    )
                    Computation._repeat_output(
                        output,
                        Computation._output_paths_for_each_repetition_count(
                            output_paths_with_repetition_counts
                        ),
                    )
                else:
                    assert len(output_paths_with_repetition_counts) == 1, (
                        "The set of output paths with repetition counts of the function '"
                        + str(self)
                        + "', namely '"
                        + str(output_paths_with_repetition_counts)
                        + "' contains more than 1 element"
                    )
                    # `next(iter(set))` returns one element of a set without removing it
                    output_path, _ = next(iter(output_paths_with_repetition_counts))
                    parameter_name_to_value = {
                        **input_parameter_name_to_path,
                        **{"output_path": output_path},
                    }
                    output = self.function(**parameter_name_to_value)
                    assert output is None, (
                        "The function '"
                        + str(self)
                        + "' returned output although it is supposed to write output to the output path '"
                        + str(output_path)
                        + "' given as parameter 'output_path'"
                    )

        return kallable


class SaveToFile(Computation):
    """Representation of a computation that writes inputs to files."""

    def __init__(
        self, label_to_input_with_path: Dict[DataFlowLabel, Tuple[Input, Path]]
    ) -> None:
        """The constructor initializes a computation that writes inputs to files.

        Args:
            label_to_input_with_path: The map from labels to tuples of input
                and path to which the input is to be written. If the input is
                not a series, then it's gonna be written to a file with the
                given path; and, if the input is a series, then it's gonna
                be written to files within the given path, which is gonna be
                a directory. Note that series input must have ``'*'`` as
                wildcard.

        Raises:
            AssertionError: If at least one series input does not have ``'*'`` as
                wildcard.

        Examples:
            >>> save_file = SaveToFile({
            ...     DataFlowLabel("my"): (Input(data_format=DataFormat.TXT), Path("./my.txt")),
            ...     DataFlowLabel("yours"): (Input(data_format=DataFormat.VEC), Path("./yours.vec")),
            ... })

            >>> save_series = SaveToFile({
            ...     DataFlowLabel("ours"): (
            ...         Input(
            ...             data_format=DataFormat.HDR,
            ...             transmission_kind=TransmissionKind.FILE,
            ...             series=Series(format="%03d", wildcard="*")
            ...         ),
            ...         Path("./my/"),
            ...     ),
            ... })

            >>> SaveToFile({
            ...     DataFlowLabel("ours"): (
            ...         Input(
            ...             data_format=DataFormat.HDR,
            ...             transmission_kind=TransmissionKind.FILE,
            ...             series=Series(format="%03d")
            ...         ),
            ...         Path("./my/"),
            ...     ),
            ... })
            Traceback (most recent call last):
            AssertionError: The wildcard '%03d' of the series input 'Input(data_format: DataFormat.HDR, series: Series(format=%03d, wildcard=%03d), transmission_kind: TransmissionKind.FILE, repetition_count: 1)' with the output path 'my' is not '*'.
        """
        label_to_input = {}
        label_to_path = {}
        for label, input_with_path in label_to_input_with_path.items():
            inpud, path = input_with_path
            if inpud.series is not None:
                assert inpud.series.wildcard == "*", (
                    "The wildcard '"
                    + str(inpud.series.wildcard)
                    + "' of the series input '"
                    + str(inpud)
                    + "' with the output path '"
                    + str(path)
                    + "' is not '*'."
                )
            label_to_input[label] = inpud
            label_to_path[label] = path
        super().__init__(label_to_input=label_to_input, output=None)
        self.label_to_path = label_to_path

    def _members(self) -> Tuple[Any, ...]:
        """The members that uniquely identify instances of ``self``."""
        return super()._members() + (self.label_to_path,)

    def __str__(self) -> str:
        return "SaveToFile(label_to_path={label_to_path}) < {computation}".format(
            label_to_path=self.label_to_path, computation=super().__str__()
        )

    @property
    def meta_data(self) -> Dict[str, Any]:
        """See ``Computation.meta_data``"""
        return {**super().meta_data, **{"label_to_path": self.label_to_path}}

    def with_file_inputs_for(self, input_labels: Set[DataFlowLabel]) -> "SaveToFile":
        """Returns itself"""
        return self

    def create_callable(
        self, computation_label: ComputationLabel, *, debug: bool = False
    ) -> Callable[
        [
            NamedArg(Dict[DataFlowLabel, Path], "label_to_input_path"),
            NamedArg(Set[Tuple[Path, int]], "output_paths_with_repetition_counts"),
        ],
        None,
    ]:
        """See ``Computation.create_callable``"""

        def kallable(  # pylint: disable=invalid-name
            *,
            label_to_input_path: Dict[DataFlowLabel, Path],
            output_paths_with_repetition_counts: Set[Tuple[Path, int]]
        ) -> None:
            self._assert_existence_of_output_paths(output_paths_with_repetition_counts)
            print("Saving inputs of {} to files".format(computation_label))
            for input_label, input_path in label_to_input_path.items():
                output_path = self.label_to_path[input_label]
                inpud = self.input(input_label)
                if inpud.series is None:
                    if output_path.exists():
                        output_path.unlink()
                    else:
                        output_path.parent.mkdir(parents=True, exist_ok=True)
                else:
                    if output_path.exists():
                        for file_path in output_path.glob("*"):
                            file_path.unlink()
                    else:
                        output_path.mkdir(parents=True)
                Computation.error_reporting_run(
                    "cp {inpud} {output}".format(inpud=input_path, output=output_path)
                )

        return kallable


class Helper:
    @staticmethod
    def determine_repetition_counts(
        computation: Computation,
        *,
        label_to_file_input_path: Dict[DataFlowLabel, Path],
        working_directory_path: Optional[Path] = None,
        debug: bool = False
    ) -> Dict[DataFlowLabel, int]:
        """Tries to determine the repetition counts of all inputs of the given
        computation.

        Note that some commands consume inputs multiple times although they are
        only needed once. For example, the command ``rfluxmtx -v -I+ -ab 4 -ad
        5000 -lw 0.0002 -n 16 -y 100 - {geometry} -i {octree} < {points}``
        consumes the input ``points`` as often as it is supplied although
        supplying it only once and waiting for the command to finish works just
        fine. Therefore, you should do at least a sanity check of the
        repetition counts determined by this function and maybe even some
        custom tests to check whether the counts are actually correct.

        Args:
            computation: The computation for which the input repetition counts
                shall be determined.
            label_to_file_input_path: A map from input labels to file paths or
                file series paths with wildcards according to the expected
                input meta data of the computation.
            working_directory_path: The path to the directory in which
                intermediary results are stored. If given, you must make
                sure that the directory is empty. If ``None``, an empty temporary
                directory is created. If ``debug`` is ``False``, once execution
                finishes, the working directory and its content is removed.
            debug: Whether to print debug output and don't clean up the working
                directory. Is ``False`` by default.

        Returns:
            A map from input labels to the determined repetition counts that
            are hopefully correct.

        Raises:
            AssertionError: If for a file input of the computation, the given
                path does not designate a file.
            AssertionError: If for a series input of the computation, the given
                path's directory does not designate a directory or the path
                does not contain the series' wildcard.

        Examples:
            >>> tmp = getfixture('tmp_path')

            >>> txt_file_path = tmp / "my.txt"
            >>> txt_file_path.write_text("1 1 2 3 5 8 13 21 ...")
            21

            >>> count_words = ShellCommand(
            ...     "wc --words {a}",
            ...     label_to_input={DataFlowLabel("a"): Input(data_format=DataFormat.TXT),},
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> label_to_count = Helper.determine_repetition_counts(
            ...     count_words,
            ...     label_to_file_input_path={DataFlowLabel("a"): txt_file_path},
            ... )
            >>> {str(label): count for label, count in label_to_count.items()}
            {'a': 1}

            >>> count_bytes_and_chars_and_words = ShellCommand(
            ...     "wc --bytes {a} && wc --chars {a} && wc --words {a}",
            ...     label_to_input={DataFlowLabel("a"): Input(data_format=DataFormat.TXT),},
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> label_to_count = Helper.determine_repetition_counts(
            ...     count_bytes_and_chars_and_words,
            ...     label_to_file_input_path={DataFlowLabel("a"): txt_file_path},
            ... )
            >>> {str(label): count for label, count in label_to_count.items()}
            {'a': 3}

            >>> another_txt_file_path = tmp / "my.txt"
            >>> another_txt_file_path.write_text("2 3 5 7 11 13 17 19 23 ...")
            26

            >>> count_lines_and_merge_two_lists = ShellCommand(
            ...     "wc --lines {first_list} && paste {first_list} {second_list}",
            ...     label_to_input={
            ...         DataFlowLabel("first_list"): Input(data_format=DataFormat.TXT),
            ...         DataFlowLabel("second_list"): Input(data_format=DataFormat.TXT),
            ...     },
            ...     output=Output(data_format=DataFormat.TXT),
            ... )
            >>> label_to_count = Helper.determine_repetition_counts(
            ...     count_lines_and_merge_two_lists,
            ...     label_to_file_input_path={
            ...         DataFlowLabel("first_list"): txt_file_path,
            ...         DataFlowLabel("second_list"): another_txt_file_path,
            ...     },
            ... )
            >>> {str(label): count for label, count in label_to_count.items()}
            {'first_list': 2, 'second_list': 1}

            >>> series_directory_path = tmp / "series"
            >>> series_directory_path.mkdir()
            >>> (series_directory_path / "list-1.txt").write_text("1 2 3 4")
            7
            >>> (series_directory_path / "list-2.txt").write_text("2 4 6 8")
            7
            >>> (series_directory_path / "list-3.txt").write_text("4 8 12 16")
            9

            >>> merge_multiple_lists = ShellCommand(
            ...     "wc --words {zeroeth_list} && wc --lines {first_list} && paste {multiple_lists} && touch {output}",
            ...     label_to_input={
            ...         DataFlowLabel("zeroeth_list"): Input(
            ...             data_format=DataFormat.TXT,
            ...             transmission_kind=TransmissionKind.FILE,
            ...         ),
            ...         DataFlowLabel("first_list"): Input(data_format=DataFormat.TXT),
            ...         DataFlowLabel("multiple_lists"): Input(
            ...             data_format=DataFormat.TXT,
            ...             transmission_kind=TransmissionKind.FILE,
            ...             series=Series(format="%03d", wildcard="*"),
            ...         ),
            ...     },
            ...     output=Output(data_format=DataFormat.TXT, series=Series(format="%03d")),
            ... )
            >>> label_to_count = Helper.determine_repetition_counts(
            ...     merge_multiple_lists,
            ...     label_to_file_input_path={
            ...         DataFlowLabel("zeroeth_list"): txt_file_path,
            ...         DataFlowLabel("first_list"): another_txt_file_path,
            ...         DataFlowLabel("multiple_lists"): series_directory_path / "list-*.txt",
            ...     },
            ... )
            >>> {str(label): count for label, count in label_to_count.items()}
            {'zeroeth_list': 1, 'first_list': 1, 'multiple_lists': 1}

            >>> Helper.determine_repetition_counts(
            ...     count_lines_and_merge_two_lists,
            ...     label_to_file_input_path={
            ...         DataFlowLabel("first_list"): tmp / "non-existent.file",
            ...         DataFlowLabel("second_list"): another_txt_file_path,
            ...     },
            ... ) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: There is nothing under the input path '.../non-existent.file'.

            >>> directory_path = tmp / "directory"
            >>> directory_path.mkdir()
            >>> Helper.determine_repetition_counts(
            ...     count_lines_and_merge_two_lists,
            ...     label_to_file_input_path={
            ...         DataFlowLabel("first_list"): directory_path,
            ...         DataFlowLabel("second_list"): another_txt_file_path,
            ...     },
            ... ) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: Whatever the input path '.../directory' designates is not a file.

            >>> Helper.determine_repetition_counts(
            ...     merge_multiple_lists,
            ...     label_to_file_input_path={
            ...         DataFlowLabel("zeroeth_list"): txt_file_path,
            ...         DataFlowLabel("first_list"): another_txt_file_path,
            ...         DataFlowLabel("multiple_lists"): series_directory_path / "list-%03d.txt",
            ...     },
            ... ) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The wildcard '*' does not occur in the input path '.../series/list-%03d.txt'.
        """
        try:
            actual_working_directory_path = (
                Path(mkdtemp())
                if working_directory_path is None
                else working_directory_path
            )
            label_to_input_path = Helper._label_to_input_path(
                computation=computation,
                label_to_file_input_path=label_to_file_input_path,
                working_directory_path=actual_working_directory_path,
            )
            output_and_path = Helper._setup_output_path(
                computation=computation,
                working_directory_path=actual_working_directory_path,
            )
            process = Process(
                target=computation.create_callable(
                    ComputationLabel(("computation",)), debug=debug
                ),
                kwargs={
                    "label_to_input_path": label_to_input_path,
                    "output_paths_with_repetition_counts": set()
                    if output_and_path is None
                    else {(output_and_path[1], 1)},
                },
            )
            process.start()
            return Helper._pipe_inputs_repeatedly_until_process_finished(
                computation=computation,
                process=process,
                label_to_file_input_path=label_to_file_input_path,
                label_to_input_path=label_to_input_path,
            )
        finally:
            if not debug and actual_working_directory_path.exists():
                shutil.rmtree(str(actual_working_directory_path))

    @staticmethod
    def _label_to_input_path(
        computation: Computation,
        label_to_file_input_path: Dict[DataFlowLabel, Path],
        working_directory_path: Path,
    ) -> Dict[DataFlowLabel, Path]:
        result = {}
        for label, input_path in label_to_file_input_path.items():
            inpud = computation.input(label)
            if inpud.series is None:
                assert input_path.exists(), (
                    "There is nothing under the input path '" + str(input_path) + "'."
                )
                assert input_path.is_file(), (
                    "Whatever the input path '"
                    + str(input_path)
                    + "' designates is not a file."
                )
            else:
                assert input_path.parent.exists(), (
                    "There is nothing under the path '" + str(input_path.parent) + "'."
                )
                assert input_path.parent.is_dir(), (
                    "Whatever the path '"
                    + str(input_path.parent)
                    + "' designates is not a directory."
                )
                assert inpud.series.wildcard in str(input_path), (
                    "The wildcard '"
                    + str(inpud.series.wildcard)
                    + "' does not occur in the input path '"
                    + str(input_path)
                    + "'."
                )
            if inpud.transmission_kind == TransmissionKind.FILE:
                result[label] = input_path
            elif inpud.transmission_kind == TransmissionKind.NAMED_PIPE:
                path = Path(
                    working_directory_path,
                    "pipe_for_{label}.{extension}".format(
                        label=label, extension=inpud.file_extension()
                    ),
                )
                os.mkfifo(str(path))
                result[label] = path
            else:
                raise NotImplementedError
        return result

    @staticmethod
    def _setup_output_path(
        computation: Computation, working_directory_path: Path
    ) -> Optional[Tuple[Output, Path]]:
        output = computation.output
        if output is None:
            return None
        if output.series is None:
            output_path = working_directory_path / "output.{extension}".format(
                extension=output.file_extension()
            )
            output_path.touch()
        else:
            output_path = (
                working_directory_path
                / "output"
                / "file_{wildcard}.{extension}".format(
                    wildcard=output.series.wildcard, extension=output.file_extension()
                )
            )
            output_path.parent.mkdir()
        return (output, output_path)

    @staticmethod
    def _pipe_inputs_repeatedly_until_process_finished(
        computation: Computation,
        process: Process,
        label_to_file_input_path: Dict[DataFlowLabel, Path],
        label_to_input_path: Dict[DataFlowLabel, Path],
    ) -> Dict[DataFlowLabel, int]:
        label_to_process = {}  # type: Dict[DataFlowLabel, Process]
        label_to_repetition_count = {
            label: (
                0
                if computation.input(label).transmission_kind
                == TransmissionKind.NAMED_PIPE
                else 1
            )
            for label in label_to_input_path
        }
        while process.is_alive():
            for label, file_input_path in label_to_file_input_path.items():
                inpud = computation.input(label)
                if inpud.transmission_kind == TransmissionKind.NAMED_PIPE:
                    if (
                        label not in label_to_process
                        or not label_to_process[label].is_alive()
                    ):
                        input_process = Process(
                            target=Computation.error_reporting_run,
                            args=(
                                "cat {inpud} > {output}".format(
                                    inpud=file_input_path,
                                    output=label_to_input_path[label],
                                ),
                            ),
                        )
                        input_process.start()
                        label_to_process[label] = input_process
                        label_to_repetition_count[label] += 1

        process.join()

        for label, input_process in label_to_process.items():
            if input_process.is_alive():
                label_to_repetition_count[label] -= 1
                input_process.terminate()  # TODO In Python 3.7 use `close` (https://docs.python.org/3/library/multiprocessing.html#multiprocessing.Process.close)

        return label_to_repetition_count


"""Data flows"""  # pylint: disable=pointless-string-statement
DataFlow = Tuple[Computation, Computation]

"""Data-flow graphs"""  # pylint: disable=pointless-string-statement
DataFlowGraph = LabeledDirectedAcyclicGraph[
    Computation, ComputationLabel, DataFlowLabel
]
