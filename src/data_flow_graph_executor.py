"""Data-flow graph executors execute the computations defined by a data-flow
graph in parallel whenever possible, where data flows from child to parent
computations through named pipes whenever possible.
"""

import time
import shutil
from labeled_directed_acyclic_graph import LabeledDirectedAcyclicGraph
from data_flow_graph import (
    Computation,
    ComputationLabel,
    DataFlow,
    DataFlowGraph,
    DataFlowLabel,
    DataFormat,
    File,
    File,
    Function,
    Input,
    InputOrOutput,
    Output,
    SaveToFile,
    Series,
    ShellCommand,
    TransmissionKind,
)
from enum import Enum, unique
from collections import Counter
from multiprocessing import Process
import os
from tempfile import mkdtemp
from typing import Any, Dict, Iterable, KeysView, Optional, Set, Tuple, cast
from pathlib import Path
from graphviz import Digraph


@unique
class _ProcessStatus(Enum):
    """Enumeration of possible process states other than 'finished'.

    Examples:
        >>> _ProcessStatus.PENDING
        <_ProcessStatus.PENDING: 'pending'>
    """

    PENDING = "pending"
    STARTED = "started"


@unique
class Transformation(Enum):
    """Enumeration of possible transformations to apply to the data-flow graph
    before its execution to avoid deadlocks auto-magically.

    * ``Transformation.NONE`` stands for __no__ transformation.
    * ``Transformation.DIAMONDS`` stands for the transformation in which each tip
      of a diamond-shaped sub-graph is made to transmit data via files. In a
      diamond-shaped sub-graph, data flows from a source, the tip, along two
      distinct paths to a sink, another computation.
    * ``Transformation.SPLITS`` stands for the transformation in which each
      computation that transmits data to at least two other computations is
      made to transmit data via files. This includes all tips of diamond-shaped
      sub-graphs.
    * ``Transformation.ALL`` stands for the transformation in which each
      computation is made to transmit data via files. This includes all
      "splits".

    For further information see the documentation of ``DataFlowGraphExecutor``.

    Examples:
        >>> Transformation.DIAMONDS
        <Transformation.DIAMONDS: 'diamonds'>
    """

    NONE = "none"
    DIAMONDS = "diamonds"
    SPLITS = "splits"
    ALL = "all"


class DataFlowGraphExecutor:
    """Data-flow graph executors execute the computations defined by
    a data-flow graph in parallel whenever possible, where data flows from
    child to parent computations through named pipes whenever possible.

    Each computation of the data-flow graph runs in its own process and is
    started as soon as possible. If all inputs of a computation are streamed,
    then it is started immediately after its child computations (recall that
    data flows from children to parents). If some inputs are transfered as
    files, then the computation is started as soon as the respective child
    computations have finished and thus the files contain the required data.

    To maximize throughput as much data as possible should be streamed.
    Streaming is sometimes not possible because, for example, some computation
    may read input an in advance unknown number of times and other computations
    sharing some inputs may block each others input streams by reading inputs
    in different orders -- such a blockage is known as 'deadlock'.

    You can follow the execution progress by opening the visual representation
    ``data_flow_graph.gv.svg`` of the data-flow graph, for example, in a web
    browser. Pending, that is, not yet started, computations are blue, running
    computations red, and finished computations green; note that computations
    may be pending because they need child computations from which they
    receive inputs as files to finish first. Roll over a computation to get its
    details in a pop-up box and roll over a data-flow arrow to get its details.

    Examples:
        To warm up, let's execute the following self-contained shell command:

        >>> tmp = getfixture('tmp_path')
        >>> def hello_x_graph(*, x: str, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("hello_x",)),
        ...         ShellCommand(
        ...             "echo Hello {x}! > {output_path}".format(x=x, output_path=output_path),
        ...             label_to_input={},
        ...             output=None,
        ...         ),
        ...         {},
        ...     )
        >>> hello_x_output_path = tmp / "hello_x.txt"
        >>> DataFlowGraphExecutor(hello_x_graph(x="Bruno", output_path=hello_x_output_path)).execute()
        >>> hello_x_output_path.read_text()
        'Hello Bruno!\\n'

        Instead of hard-coding ``x`` and ``output_path`` in the shell command,
        let's read the former from a file and write the output to ``stdout``.

        >>> def y_graph(*, path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("y",)),
        ...         File(path),
        ...         {},
        ...     )
        >>> def hello_y_graph(*, y_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("hello_y",)),
        ...         ShellCommand(
        ...             'echo Hello "$(cat {y})"!',
        ...             label_to_input={DataFlowLabel("y"): Input(data_format=DataFormat.TXT)},
        ...             output=Output(data_format=DataFormat.TXT),
        ...         ),
        ...         {DataFlowLabel("y"): y_graph(path=y_path)},
        ...     )
        >>> def save_hello_y_graph(*, y_path: Path, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("save_hello_y_graph",)),
        ...         SaveToFile({
        ...             DataFlowLabel("hello_y"): (Input(data_format=DataFormat.TXT), output_path),
        ...         }),
        ...         {DataFlowLabel("hello_y"): hello_y_graph(y_path=y_path)},
        ...     )
        >>> y_path = tmp / "y.txt"
        >>> y_path.write_text("Simon")
        5
        >>> hello_y_output_path = tmp / "hello_y.txt"
        >>> DataFlowGraphExecutor(save_hello_y_graph(y_path=y_path, output_path=hello_y_output_path)).execute()
        >>> hello_y_output_path.read_text()
        'Hello Simon!\\n'

        Before saving the output, let's add an additional transformation step.

        >>> def goodbye_y_graph(*, y_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("goodbye_y",)),
        ...         ShellCommand(
        ...             'sed --expression "s/Hello/Goodbye/" {hello_y}',
        ...             label_to_input={DataFlowLabel("hello_y"): Input(data_format=DataFormat.TXT)},
        ...             output=Output(data_format=DataFormat.TXT),
        ...         ),
        ...         {DataFlowLabel("hello_y"): hello_y_graph(y_path=y_path)},
        ...     )
        >>> def save_goodbye_y_graph(*, y_path: Path, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("save_goodbye_y_graph",)),
        ...         SaveToFile({
        ...             DataFlowLabel("goodbye_y"): (Input(data_format=DataFormat.TXT), output_path),
        ...         }),
        ...         {DataFlowLabel("goodbye_y"): goodbye_y_graph(y_path=y_path)},
        ...     )
        >>> y_path.write_text("Johannes")
        8
        >>> goodbye_y_output_path = tmp / "goodbye_y.txt"
        >>> DataFlowGraphExecutor(save_goodbye_y_graph(y_path=y_path, output_path=goodbye_y_output_path)).execute()
        >>> goodbye_y_output_path.read_text()
        'Goodbye Johannes!\\n'

        For a change, let's do a transformation with a good old Python function.

        >>> def see_ya_y_graph(*, y_path: Path) -> DataFlowGraph:
        ...     def see_ya_y(hello_y_input_path: Path) -> Optional[bytes]:
        ...         return hello_y_input_path.read_text().replace("Hello", "See ya,").encode("utf-8")
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("see_ya_y",)),
        ...         Function(
        ...             see_ya_y,
        ...             label_to_input={DataFlowLabel("hello_y"): Input(data_format=DataFormat.TXT)},
        ...             output=Output(data_format=DataFormat.TXT),
        ...         ),
        ...         {DataFlowLabel("hello_y"): hello_y_graph(y_path=y_path)},
        ...     )
        >>> def save_see_ya_y_graph(*, y_path: Path, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("save_see_ya_y",)),
        ...         SaveToFile({
        ...             DataFlowLabel("see_ya_y"): (Input(data_format=DataFormat.TXT), output_path),
        ...         }),
        ...         {DataFlowLabel("see_ya_y"): see_ya_y_graph(y_path=y_path)},
        ...     )
        >>> y_path.write_text("Jan")
        3
        >>> see_ya_y_output_path = tmp / "see_ya_y.txt"
        >>> DataFlowGraphExecutor(save_see_ya_y_graph(y_path=y_path, output_path=see_ya_y_output_path)).execute()
        >>> see_ya_y_output_path.read_text()
        'See ya, Jan!\\n'

        Now, let's read from one computation twice and combine the outputs of
        two computations (have a look at the data-flow graph after executing
        it by opening ``data_flow_graph.gv.svg`` in a web browser).

        >>> def hello_and_see_ya_y_graph(*, y_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("hello_and_see_ya_y",)),
        ...         ShellCommand(
        ...             "cat {hello_y} {hello_y} {see_ya_y}",
        ...             label_to_input={
        ...                 DataFlowLabel("hello_y"): Input(data_format=DataFormat.TXT, repetition_count=2),
        ...                 DataFlowLabel("see_ya_y"): Input(data_format=DataFormat.TXT),
        ...             },
        ...             output=Output(data_format=DataFormat.TXT),
        ...         ),
        ...         {
        ...             DataFlowLabel("hello_y"): hello_y_graph(y_path=y_path),
        ...             DataFlowLabel("see_ya_y"): see_ya_y_graph(y_path=y_path),
        ...         },
        ...     )
        >>> def save_hello_and_see_ya_y_graph(*, y_path: Path, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("save_hello_and_see_ya_y",)),
        ...         SaveToFile({
        ...             DataFlowLabel("hello_and_see_ya_y"): (Input(data_format=DataFormat.TXT), output_path),
        ...         }),
        ...         {DataFlowLabel("hello_and_see_ya_y"): hello_and_see_ya_y_graph(y_path=y_path)},
        ...     )
        >>> y_path.write_text("Adrienne")
        8
        >>> hello_and_see_ya_y_output_path = tmp / "hello_and_see_ya_y.txt"
        >>> DataFlowGraphExecutor(save_hello_and_see_ya_y_graph(y_path=y_path, output_path=hello_and_see_ya_y_output_path)).execute()
        >>> hello_and_see_ya_y_output_path.read_text()
        'Hello Adrienne!\\nHello Adrienne!\\nSee ya, Adrienne!\\n'

        Finally, let's process a file series (this is what the Radiance
        commands ``rfluxmtx`` and ``dctimestep`` and others consume and produce).

        >>> def names_graph(*, path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("names",)),
        ...         File(path),
        ...         {},
        ...     )
        >>> def wish_good_mornings_graph(*, names_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("wish_good_mornings",)),
        ...         ShellCommand(
        ...             formatted_command=(
        ...                 "for idx in {{0..9}}\\n"
        ...                 "do\\n"
        ...                 "    if test -f \\"{names}\\"; then\\n"
        ...                 "       echo \\"Good morning $(cat {names})!\\" > {output}\\n"
        ...                 "    fi\\n"
        ...                 "done"
        ...             ),
        ...             label_to_input={
        ...                 DataFlowLabel("names"): Input(data_format=DataFormat.TXT, transmission_kind=TransmissionKind.FILE, series=Series(format="%01d", wildcard="${idx}"))
        ...             },
        ...             output=Output(data_format=DataFormat.TXT, series=Series(format="%01d", wildcard="${idx}")),
        ...         ),
        ...         {DataFlowLabel("names"): names_graph(path=names_path)},
        ...     )
        >>> def wish_good_evenings_graph(*, names_path: Path) -> DataFlowGraph:
        ...     def wish_good_evenings(names_input_path: Path, output_path: Path) -> Optional[bytes]:
        ...         for idx in range(9):
        ...             name_input_path = Path(str(names_input_path).format(idx=idx))
        ...             if name_input_path.exists():
        ...                 name = name_input_path.read_text()
        ...                 Path(str(output_path).format(idx=idx)).write_text("Good evening {}!".format(name))
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("wish_good_evenings",)),
        ...         Function(
        ...             wish_good_evenings,
        ...             label_to_input={
        ...                 DataFlowLabel("names"): Input(data_format=DataFormat.TXT, transmission_kind=TransmissionKind.FILE, series=Series(format="%01d", wildcard="{idx}"))
        ...             },
        ...             output=Output(data_format=DataFormat.TXT, series=Series(format="%01d", wildcard="{idx}")),
        ...         ),
        ...         {DataFlowLabel("names"): names_graph(path=names_path)},
        ...     )
        >>> def wish_goods_graph(*, names_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("wish_goods",)),
        ...         ShellCommand(
        ...             formatted_command=(
        ...                 "for idx in {{0..9}}\\n"
        ...                 "do\\n"
        ...                 "    if test -f \\"{mornings}\\" && test -f \\"{evenings}\\"; then\\n"
        ...                 "       cat {mornings} {evenings} > {output}\\n"
        ...                 "    fi\\n"
        ...                 "done"
        ...             ),
        ...             label_to_input={
        ...                 DataFlowLabel("mornings"): Input(data_format=DataFormat.TXT, transmission_kind=TransmissionKind.FILE, series=Series(format="%01d", wildcard="${idx}")),
        ...                 DataFlowLabel("evenings"): Input(data_format=DataFormat.TXT, transmission_kind=TransmissionKind.FILE, series=Series(format="%01d", wildcard="${idx}")),
        ...             },
        ...             output=Output(data_format=DataFormat.TXT, series=Series(format="%01d", wildcard="${idx}")),
        ...         ),
        ...         {
        ...             DataFlowLabel("mornings"): wish_good_mornings_graph(names_path=names_path),
        ...             DataFlowLabel("evenings"): wish_good_evenings_graph(names_path=names_path),
        ...         },
        ...     )
        >>> def save_wish_goods_graph(*, names_path: Path, output_path: Path) -> DataFlowGraph:
        ...     return DataFlowGraph.construct(
        ...         ComputationLabel(("save_wish_goods",)),
        ...         SaveToFile({
        ...             DataFlowLabel("wish_goods"): (Input(data_format=DataFormat.TXT, transmission_kind=TransmissionKind.FILE, series=Series(format="%01d", wildcard="*")), output_path),
        ...         }),
        ...         {DataFlowLabel("wish_goods"): wish_goods_graph(names_path=names_path)},
        ...     )
        >>> names = ["Bruno", "Simon", "Johannes", "Jan", "Adriennce", "Christoph"]
        >>> names_directory_path = tmp / "names"
        >>> names_directory_path.mkdir()
        >>> for i, name in enumerate(names):
        ...     _ = (names_directory_path / "name-{}.txt".format(i)).write_text(name)
        >>> names_path = names_directory_path / "name-%01d.txt"
        >>> wish_goods_output_path = tmp / "wish_goods"
        >>> DataFlowGraphExecutor(save_wish_goods_graph(names_path=names_path, output_path=wish_goods_output_path)).execute()
        >>> wish_goods_output_path.is_dir()
        True
        >>> {str(path.name) for path in wish_goods_output_path.iterdir()} == {'file_0.txt', 'file_1.txt', 'file_2.txt', 'file_3.txt', 'file_4.txt', 'file_5.txt'}
        True
        >>> [(wish_goods_output_path / "file_{}.txt".format(i)).read_text() for i, _ in enumerate(names)]
        ['Good morning Bruno!\\nGood evening Bruno!', 'Good morning Simon!\\nGood evening Simon!', 'Good morning Johannes!\\nGood evening Johannes!', 'Good morning Jan!\\nGood evening Jan!', 'Good morning Adriennce!\\nGood evening Adriennce!', 'Good morning Christoph!\\nGood evening Christoph!']
    """

    @staticmethod
    def _transform_the_data_flow_graph_pessimistically_to_avoid_deadlocks(
        data_flow_graph: DataFlowGraph, *, transformation: Transformation
    ) -> DataFlowGraph:
        # TODO All transformations have the same frame: Make a meta
        #      transformation with that frame that is being passed the individual
        #      transformation approach as callable; this removes duplication
        #      and makes the code easier to digest.
        if transformation == Transformation.NONE:
            return data_flow_graph
        if transformation == Transformation.DIAMONDS:
            return DataFlowGraphExecutor._transform_diamond_tips_pessimistically_to_avoid_deadlocks(
                data_flow_graph
            )
        if transformation == Transformation.SPLITS:
            return DataFlowGraphExecutor._transform_splits_pessimistically_to_avoid_deadlocks(
                data_flow_graph
            )
        if transformation == Transformation.ALL:
            return DataFlowGraphExecutor._transform_all_computations_pessimistically_to_avoid_deadlocks(
                data_flow_graph
            )
        raise NotImplementedError

    @staticmethod
    def _transform_diamond_tips_pessimistically_to_avoid_deadlocks(
        data_flow_graph: DataFlowGraph
    ) -> DataFlowGraph:
        # pylint: disable=too-many-nested-blocks
        result = LabeledDirectedAcyclicGraph()  # type: DataFlowGraph
        computation_to_inputs_that_need_to_be_files = (  # pylint: disable=invalid-name
            {}
        )  # type: Dict[Computation, Set[DataFlowLabel]]
        for computation in data_flow_graph.traverse_nodes_layerwise_deepest_first():
            possibly_transformed_computation = computation.with_file_inputs_for(  # pylint: disable=invalid-name
                computation_to_inputs_that_need_to_be_files.get(computation, set())
            )
            result.add_node(
                data_flow_graph.node_label(computation),
                possibly_transformed_computation,
                data_flow_graph.edge_label_to_child_node_label(computation),
            )
            if not isinstance(computation, File):
                pipe_input_parent_computations = [
                    parent_computation
                    for parent_computation in data_flow_graph.parent_nodes(computation)
                    if parent_computation.input(
                        data_flow_graph.edge_label((parent_computation, computation))
                    ).transmission_kind
                    == TransmissionKind.NAMED_PIPE
                ]
                if len(pipe_input_parent_computations) >= 2:
                    ancestors = [
                        data_flow_graph.ancestors(parent_computation)
                        for parent_computation in pipe_input_parent_computations
                    ]
                    for i, ith_computation in enumerate(pipe_input_parent_computations):
                        for k, kth_computation in enumerate(
                            pipe_input_parent_computations[i + 1 :]
                        ):
                            if not ancestors[i].isdisjoint(ancestors[k]):
                                if (
                                    ith_computation
                                    not in computation_to_inputs_that_need_to_be_files
                                ):
                                    computation_to_inputs_that_need_to_be_files[
                                        ith_computation
                                    ] = set()
                                computation_to_inputs_that_need_to_be_files[
                                    ith_computation
                                ].add(
                                    data_flow_graph.edge_label(
                                        (ith_computation, computation)
                                    )
                                )
                                if (
                                    kth_computation
                                    not in computation_to_inputs_that_need_to_be_files
                                ):
                                    computation_to_inputs_that_need_to_be_files[
                                        kth_computation
                                    ] = set()
                                computation_to_inputs_that_need_to_be_files[
                                    kth_computation
                                ].add(
                                    data_flow_graph.edge_label(
                                        (kth_computation, computation)
                                    )
                                )
        return result

    @staticmethod
    def _transform_splits_pessimistically_to_avoid_deadlocks(
        data_flow_graph: DataFlowGraph
    ) -> DataFlowGraph:
        # pylint: disable=too-many-nested-blocks
        result = LabeledDirectedAcyclicGraph()  # type: DataFlowGraph
        computation_to_inputs_that_need_to_be_files = (  # pylint: disable=invalid-name
            {}
        )  # type: Dict[Computation, Set[DataFlowLabel]]
        for computation in data_flow_graph.traverse_nodes_layerwise_deepest_first():
            possibly_transformed_computation = computation.with_file_inputs_for(  # pylint: disable=invalid-name
                computation_to_inputs_that_need_to_be_files.get(computation, set())
            )
            result.add_node(
                data_flow_graph.node_label(computation),
                possibly_transformed_computation,
                data_flow_graph.edge_label_to_child_node_label(computation),
            )
            if not isinstance(computation, File):
                pipe_input_parent_computations = [
                    parent_computation
                    for parent_computation in data_flow_graph.parent_nodes(computation)
                    if parent_computation.input(
                        data_flow_graph.edge_label((parent_computation, computation))
                    ).transmission_kind
                    == TransmissionKind.NAMED_PIPE
                ]
                if len(pipe_input_parent_computations) >= 2:
                    for parent_computation in pipe_input_parent_computations:
                        if (
                            parent_computation
                            not in computation_to_inputs_that_need_to_be_files
                        ):
                            computation_to_inputs_that_need_to_be_files[
                                parent_computation
                            ] = set()
                        computation_to_inputs_that_need_to_be_files[
                            parent_computation
                        ].add(
                            data_flow_graph.edge_label(
                                (parent_computation, computation)
                            )
                        )
        return result

    @staticmethod
    def _transform_all_computations_pessimistically_to_avoid_deadlocks(
        data_flow_graph: DataFlowGraph
    ) -> DataFlowGraph:
        # pylint: disable=too-many-nested-blocks
        result = LabeledDirectedAcyclicGraph()  # type: DataFlowGraph
        computation_to_inputs_that_need_to_be_files = (  # pylint: disable=invalid-name
            {}
        )  # type: Dict[Computation, Set[DataFlowLabel]]
        for computation in data_flow_graph.traverse_nodes_layerwise_deepest_first():
            possibly_transformed_computation = computation.with_file_inputs_for(  # pylint: disable=invalid-name
                computation_to_inputs_that_need_to_be_files.get(computation, set())
            )
            result.add_node(
                data_flow_graph.node_label(computation),
                possibly_transformed_computation,
                data_flow_graph.edge_label_to_child_node_label(computation),
            )
            if not isinstance(computation, File):
                pipe_input_parent_computations = [
                    parent_computation
                    for parent_computation in data_flow_graph.parent_nodes(computation)
                    if parent_computation.input(
                        data_flow_graph.edge_label((parent_computation, computation))
                    ).transmission_kind
                    == TransmissionKind.NAMED_PIPE
                ]
                for parent_computation in pipe_input_parent_computations:
                    if (
                        parent_computation
                        not in computation_to_inputs_that_need_to_be_files
                    ):
                        computation_to_inputs_that_need_to_be_files[
                            parent_computation
                        ] = set()
                    computation_to_inputs_that_need_to_be_files[parent_computation].add(
                        data_flow_graph.edge_label((parent_computation, computation))
                    )
        return result

    @staticmethod
    def _assert_validity(data_flow_graph: DataFlowGraph) -> None:
        for computation in data_flow_graph.nodes:
            DataFlowGraphExecutor._assert_that_child_computations_have_output(
                data_flow_graph, computation=computation
            )
            label_to_output = {
                data_flow_graph.edge_label((computation, child_computation)): cast(
                    Output, child_computation.output
                )
                for child_computation in data_flow_graph.child_nodes(computation)
            }
            DataFlowGraphExecutor._assert_that_data_flow_labels_and_input_labels_match(
                data_flow_graph,
                computation=computation,
                data_flow_labels=label_to_output.keys(),
            )
            DataFlowGraphExecutor._assert_that_outputs_and_inputs_match(
                data_flow_graph,
                computation=computation,
                label_to_output=label_to_output,
            )

    @staticmethod
    def _assert_that_child_computations_have_output(
        data_flow_graph: DataFlowGraph, *, computation: Computation
    ) -> None:
        for child_computation in data_flow_graph.child_nodes(computation):
            assert child_computation.output is not None, (
                "The child computation '"
                + str(data_flow_graph.node_label(child_computation))
                + "' of '"
                + str(data_flow_graph.node_label(computation))
                + "' has no output"
            )

    @staticmethod
    def _assert_that_data_flow_labels_and_input_labels_match(
        data_flow_graph: DataFlowGraph,
        *,
        computation: Computation,
        data_flow_labels: KeysView[DataFlowLabel]
    ) -> None:
        input_labels = computation.input_labels
        assert Counter(data_flow_labels) == Counter(input_labels), (
            "The data flow labels '"
            + ", ".join([str(data_flow_label) for data_flow_label in data_flow_labels])
            + "' are not the same as the input labels '"
            + ", ".join([str(input_label) for input_label in input_labels])
            + "' of the computation '"
            + str(data_flow_graph.node_label(computation))
            + "'"
        )

    @staticmethod
    def _assert_that_outputs_and_inputs_match(
        data_flow_graph: DataFlowGraph,
        *,
        computation: Computation,
        label_to_output: Dict[DataFlowLabel, Output]
    ) -> None:
        for data_flow_label, output in label_to_output.items():
            assert output.can_be_used_as_input_for(
                computation.input(data_flow_label)
            ), (
                "The actual input, namely '"
                + str(output)
                + "', for label '"
                + str(data_flow_label)
                + "' is not compatible with the expected input, namely '"
                + str(computation.input(data_flow_label))
                + "', for the same label of the computation '"
                + str(data_flow_graph.node_label(computation))
                + "'"
            )

    # Why we use named pipes instead of files: https://askubuntu.com/questions/449132/why-use-a-named-pipe-instead-of-a-file/449192#449192
    @staticmethod
    def _touch_files_and_create_named_pipes_and_create_directories(
        data_flow_graph: DataFlowGraph, *, working_directory_path: Path
    ) -> Dict[DataFlow, Path]:
        data_flow_to_path = {}
        for computation in data_flow_graph.nodes:
            for child_computation in data_flow_graph.child_nodes(computation):
                if isinstance(child_computation, File):
                    path = child_computation.path
                else:
                    path = DataFlowGraphExecutor._touch_file_or_create_named_pipe_or_create_directory(
                        data_flow_graph,
                        working_directory_path=working_directory_path,
                        computation=computation,
                        child_computation=child_computation,
                    )
                data_flow_to_path[(computation, child_computation)] = path
        return data_flow_to_path

    @staticmethod
    def _touch_file_or_create_named_pipe_or_create_directory(
        data_flow_graph: DataFlowGraph,
        *,
        working_directory_path: Path,
        computation: Computation,
        child_computation: Computation
    ) -> Path:
        computation_label = data_flow_graph.node_label(computation)
        child_label = data_flow_graph.node_label(child_computation)
        data_flow_label = data_flow_graph.edge_label((computation, child_computation))
        Path(working_directory_path, *child_label.group_path).mkdir(
            parents=True, exist_ok=True
        )
        output = cast(Output, child_computation.output)
        inpud = computation.input(data_flow_label)

        if inpud.transmission_kind == TransmissionKind.FILE:
            if inpud.series is None:
                return DataFlowGraphExecutor._touch_file(
                    working_directory_path, child_label=child_label, output=output
                )
            return DataFlowGraphExecutor._make_directory_for_series(
                working_directory_path, child_label=child_label, output=output
            )

        if inpud.transmission_kind == TransmissionKind.NAMED_PIPE:
            return DataFlowGraphExecutor._make_named_pipe(
                working_directory_path,
                child_label=child_label,
                data_flow_label=data_flow_label,
                computation_label=computation_label,
                output=output,
            )

        raise NotImplementedError

    @staticmethod
    def _touch_file(
        working_directory_path: Path, *, child_label: ComputationLabel, output: Output
    ) -> Path:
        path = Path(
            working_directory_path,
            *(child_label.group_path),
            "file_" + child_label.name + "." + output.file_extension(),
        )
        path.touch()
        return path

    @staticmethod
    def _make_named_pipe(
        working_directory_path: Path,
        *,
        child_label: ComputationLabel,
        data_flow_label: DataFlowLabel,
        computation_label: ComputationLabel,
        output: Output
    ) -> Path:
        # Note that path name follows direction of data flow not direction of edge.
        path = Path(
            working_directory_path,
            *child_label.group_path,
            "pipe_from_"
            + child_label.name
            + "_via_data_flow_"
            + str(data_flow_label)
            + "_to_"
            + computation_label.to_file_name()
            + "."
            + output.file_extension(),
        )
        os.mkfifo(str(path))
        return path

    @staticmethod
    def _make_directory_for_series(
        working_directory_path: Path, *, child_label: ComputationLabel, output: Output
    ) -> Path:
        Path(working_directory_path, *child_label).mkdir(parents=True, exist_ok=True)
        # The wildcard may differ for output and input and is thus put in when the output and input paths are created.
        return Path(
            working_directory_path,
            *child_label,
            "file_{wildcard}." + output.file_extension(),
        )

    @staticmethod
    def _fill_in_wildcard_if_necessary(
        path: Path, input_or_output: InputOrOutput
    ) -> Path:
        if input_or_output.series is None:
            return path
        return Path(str(path).format(wildcard=input_or_output.series.wildcard))

    @staticmethod
    def _terminate_alive_processes(
        processes_with_status: Iterable[Tuple[Process, _ProcessStatus]]
    ) -> None:
        for process, _ in processes_with_status:
            if process.is_alive():
                process.terminate()

    def __init__(
        self,
        data_flow_graph: DataFlowGraph,
        *,
        transformation: Transformation = Transformation.NONE
    ) -> None:
        """The constructor initializes a data-flow graph executor for the
        given data-flow graph, asserts that the graph is valid, and stores
        a copy of it with the requested transformation applied, where by
        default no transformation is applied.

        If your run into a deadlock, try the possible transformations that
        encompass each other in the following order: ``Transformation.DIAMONDS``,
        ``Transformation.SPLITS``, and ``Transformation.ALL``. From left to right,
        the possibility for deadlocks decreases to being impossible (at least
        with respect to data transmission) but the efficiency of data
        transmission decreases (because more and more transmission is made via
        files).

        Note that each executor instance is meant to be used exactly once and
        cleans up after itself properly. Preferably use it as follows:
        ``DataFlowExecutor(data_flow_graph).execute()``.

        The transformation ``Transformation.DIAMONDS`` makes deadlocks that may
        occur when data flows through named pipes in a diamond-shaped sub-graph
        impossible. Why are diamond-shaped sub-graphs problematic? Let's
        consider the case where there is a stream from a computation ``x``
        through two computations ``y1`` and ``y2`` to a computation ``z``. If ``z``
        tries to first consume the incoming data from ``y1`` completely and
        afterwards from ``y2``, (almost) nothing is going to happen at all --
        a deadlock. If ``z`` would instead consume one byte from ``y1``, then one
        from ``y2``, then one from ``y1``, and so on, then __no__ deadlock would
        occur.

        Let's see why by considering what happens in detail: First ``x`` produces
        a byte, streams it to ``y1`` and a copy to ``y2``, secondly ``y1`` and ``y2``
        consume their respective byte, start computing and say produce another
        byte each and stream it to ``z``, and thirdly ``z`` consumes the byte from
        ``y1`` but __not__ the one from ``y2``, which fills the stream buffer from
        ``y2`` to ``z`` and thus ``y2`` pauses computing waiting for the buffer to be
        drained by ``z``.

        Now, first ``x`` produces another byte, streams it to ``y1`` and a copy to
        ``y2``, secondly ``y1`` consumes its byte but __not__ ``y2`` which paused and
        thus its byte fills the stream buffer from ``x`` to ``y2`` and ``x`` pauses
        computing waiting for the buffer to be drained by ``y2``, and thirdly ``z``
        consumes the byte from ``y1`` but still __not__ the one in the buffer
        from ``y2``.

        Finally, nothing happens anymore because ``x`` cannot produce another
        byte because the stream buffer from ``x`` to ``y2`` is full, ``y2`` cannot
        consume the byte in the buffer because the one from ``y2`` to ``z`` is
        full, ``z`` waits for ``y1`` to finish before it consumes the byte in the
        buffer from ``y2`` to ``z``, and ``y1`` waits for another byte from ``x``.

        The deadlock can be avoided by having ``x`` __not__ stream data to ``y1``
        and ``y2`` byte by byte while it is computing, but instead transfer all
        bytes to ``y1`` and all bytes to ``y2`` at once through a file after it has
        finished computing. The downside is that then ``y1``, ``y2``, and ``z`` have
        to wait for ``x`` to finish before they can start computing themselves,
        and that transfer through files is less efficient than through pipes.

        Such diamond-shaped sub-graphs can be as small as in the example above
        but in general are arbitrarily big, difficult to spot, and debug
        manually. Therefore, the given data-flow graph is transformed
        automatically in a pessimistic manner to avoid such deadlocks by making
        all computations at the tip of a diamond-shape, like ``x`` above,
        transfer output as files.

        Args:
            data_flow_graph: The data-flow graph to execute.

        Raises:
            AssertionError: If the given data-flow graph is not valid, that
                is,
                * if a computation's input labels are not the same as the
                  data-flow labels to its child computations (recall that data
                  flows from children to parents), for example, if a
                  computation has ``n`` inputs and ``k`` child computations with
                  ``k != n``, or
                * if a computation has an input that is not compatible with the
                  corresponding child computation's output as defined by
                  ``Output.can_be_used_as_input_for`` (recall that data flows
                  from children to parents), for example, if the input's data
                  format is ``x`` and child's output data format is ``y`` with
                  ``y != x``.
        """
        self.debug = False
        self.working_directory_path = None  # type: Optional[Path]
        self.data_flow_to_path = {}  # type: Dict[DataFlow, Path]
        DataFlowGraphExecutor._assert_validity(data_flow_graph)
        self.data_flow_graph = DataFlowGraphExecutor._transform_the_data_flow_graph_pessimistically_to_avoid_deadlocks(
            data_flow_graph, transformation=transformation
        )

    def execute(
        self, *, working_directory_path: Optional[Path] = None, debug: bool = False
    ) -> None:
        """Executes the computations defined by a data-flow graph in parallel
        whenever possible, where data flows from child to parent computations
        through named pipes whenever possible.

        Each computation of the data-flow graph runs in its own process and is
        started as soon as possible. If all inputs of a computation are
        streamed, then it is started immediately after its child computations
        (recall that data flows from children to parents). If some inputs are
        transfered as files, then the computation is started as soon as the
        respective child computations have finished and thus the files contain
        the required data.

        To maximize throughput as much data as possible should be streamed.

        Args:
            working_directory_path: The path to the directory in which
                intermediary results are stored. If given, you must make
                sure that the directory is empty. If ``None``, an empty temporary
                directory is created. If ``debug`` is ``False``, once execution
                finishes, the working directory and its content is removed.
            debug: Whether to print debug output and don't clean up the working
                directory. Is ``False`` by default.
        """
        try:
            self.debug = debug
            self.working_directory_path = (
                Path(mkdtemp())
                if working_directory_path is None
                else working_directory_path
            )
            self.data_flow_to_path = DataFlowGraphExecutor._touch_files_and_create_named_pipes_and_create_directories(
                self.data_flow_graph, working_directory_path=self.working_directory_path
            )
            (
                computation_to_process_with_status,  # pylint: disable=invalid-name
                pending_computations,
            ) = self._init_processes()
            while (
                pending_computations
                and not self._have_all_root_computations_finished(
                    computation_to_process_with_status
                )
            ):
                started_computations = set()
                for computation in pending_computations:
                    if not self._has_to_wait_for_at_least_one_child_computation_to_finish(
                        computation, computation_to_process_with_status
                    ):
                        process, _ = computation_to_process_with_status[computation]
                        process.start()
                        computation_to_process_with_status[computation] = (
                            process,
                            _ProcessStatus.STARTED,
                        )
                        started_computations.add(computation)
                pending_computations = pending_computations.difference(
                    started_computations
                )
                _DataFlowGraphExecutionDrawer.draw(
                    self.data_flow_graph, computation_to_process_with_status
                )
                time.sleep(0.1)
            self._wait_for_root_processes_to_finish(computation_to_process_with_status)
            DataFlowGraphExecutor._terminate_alive_processes(
                computation_to_process_with_status.values()
            )
            self._assert_success_of_processes(computation_to_process_with_status)
        finally:
            self._maybe_remove_working_directory()

    def _init_processes(
        self
    ) -> Tuple[Dict[Computation, Tuple[Process, _ProcessStatus]], Set[Computation]]:
        computation_to_process_with_status = {}  # pylint: disable=invalid-name
        pending_computations = set()
        for computation in self.data_flow_graph:
            if not isinstance(computation, File):
                process = Process(
                    target=computation.create_callable(
                        self.data_flow_graph.node_label(computation), debug=self.debug
                    ),
                    kwargs={
                        "label_to_input_path": self._label_to_input_path(computation),
                        "output_paths_with_repetition_counts": self._output_paths_with_repetition_counts(
                            computation
                        ),
                    },
                )
                computation_to_process_with_status[computation] = (
                    process,
                    _ProcessStatus.PENDING,
                )
                pending_computations.add(computation)
        return (computation_to_process_with_status, pending_computations)

    def _label_to_input_path(
        self, computation: Computation
    ) -> Dict[DataFlowLabel, Path]:
        result = {}
        for child_computation in self.data_flow_graph.child_nodes(computation):
            data_flow_label = self.data_flow_graph.edge_label(
                (computation, child_computation)
            )
            result[
                data_flow_label
            ] = DataFlowGraphExecutor._fill_in_wildcard_if_necessary(
                self.data_flow_to_path[(computation, child_computation)],
                computation.input(data_flow_label),
            )
        return result

    def _output_paths_with_repetition_counts(
        self, computation: Computation
    ) -> Set[Tuple[Path, int]]:
        result = set()
        for parent_computation in self.data_flow_graph.parent_nodes(computation):
            data_flow_label = self.data_flow_graph.edge_label(
                (parent_computation, computation)
            )
            inpud = parent_computation.input(data_flow_label)
            result.add(
                (
                    DataFlowGraphExecutor._fill_in_wildcard_if_necessary(
                        self.data_flow_to_path[(parent_computation, computation)],
                        cast(Output, computation.output),
                    ),
                    inpud.repetition_count,
                )
            )
        return result

    def _have_all_root_computations_finished(
        self,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> bool:
        for (
            computation,
            process_and_status,
        ) in computation_to_process_with_status.items():
            if not self.data_flow_graph.parent_nodes(computation):
                process, status = process_and_status
                if status == _ProcessStatus.PENDING or (
                    status == _ProcessStatus.STARTED and process.is_alive()
                ):
                    return False
        return True

    def _has_to_wait_for_at_least_one_child_computation_to_finish(
        self,
        computation: Computation,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> bool:
        for child_computation in self.data_flow_graph.child_nodes(computation):
            if not isinstance(child_computation, File):
                process, status = computation_to_process_with_status[child_computation]
                data_flow_label = self.data_flow_graph.edge_label(
                    (computation, child_computation)
                )
                # Note that the returned input is not `None` due to the validity check
                if (
                    computation.input(data_flow_label).transmission_kind
                    == TransmissionKind.FILE
                ) and (
                    status == _ProcessStatus.PENDING
                    or (status == _ProcessStatus.STARTED and process.is_alive())
                ):
                    return True
        return False

    def _wait_for_root_processes_to_finish(
        self,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> None:
        # Wait for the root processes to finish (otherwise
        # `_maybe_remove_working_directory` may be called too early removing still
        # needed named pipes)
        for (
            computation,
            process_and_status,
        ) in computation_to_process_with_status.items():
            if not self.data_flow_graph.parent_nodes(computation):
                process, _ = process_and_status
                process.join()

    def _assert_success_of_processes(
        self,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> None:
        computation_labels_with_failure_exit_code = [
            (self.data_flow_graph.node_label(computation), process.exitcode)
            for computation, (
                process,
                status,
            ) in computation_to_process_with_status.items()
            if process.exitcode != 0
        ]
        assert not computation_labels_with_failure_exit_code, (
            "The following computations failed: "
            + "; ".join(
                [
                    (
                        "Computation with label '{computation_label}' exited with failure code {exit_code}."
                        " Note that the failure code `None` means that we terminated the process"
                        " because its output was not needed"
                    ).format(computation_label=computation_label, exit_code=exit_code)
                    for computation_label, exit_code in computation_labels_with_failure_exit_code
                ]
            )
        )

    def _maybe_remove_working_directory(self) -> None:
        if (
            not self.debug
            and self.working_directory_path is not None
            and self.working_directory_path.exists()
        ):
            shutil.rmtree(str(self.working_directory_path))


# ASCII graphs could be created as explained on
# https://stackoverflow.com/questions/3211801/graphviz-and-ascii-output/53895531#53895531
# Also of interest for drawing graphs could be
# http://networkx.github.io
class _DataFlowGraphExecutionDrawer:
    @staticmethod
    def draw(
        data_flow_graph: DataFlowGraph,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> None:
        # A graphical representation of the graph can be generated from a dot file with
        # `dot -Tps data_flow_graph.gv -o data_flow_graph.gvx.pdf`
        digraph = _DataFlowGraphExecutionDrawer._build_digraph()
        computation_to_label_to_input = (
            {}
        )  # type: Dict[Computation, Dict[DataFlowLabel, Input]]

        # Add computations and collect inputs
        for computation, computation_label in data_flow_graph.nodes_with_labels:
            (
                status_abbreviation,
                computation_color,
            ) = _DataFlowGraphExecutionDrawer._determine_status_abbreviation_and_color(
                computation, computation_to_process_with_status
            )
            computation_to_label_to_input[computation] = {}
            computation_input = (
                "##INPUTS##\n"
                + _DataFlowGraphExecutionDrawer._format_and_collect_inputs(
                    computation, computation_to_label_to_input[computation]
                )
            )
            computation_output = (
                "##OUTPUTS##\n"
                + _DataFlowGraphExecutionDrawer._format_output(computation.output)
            )
            meta_data = (
                "##META_DATA##\n"
                + _DataFlowGraphExecutionDrawer._format_meta_data(computation.meta_data)
            )
            tooltip = (
                computation_input + "\n\n" + computation_output + "\n\n" + meta_data
            )
            digraph.node(
                str(computation_label),
                _DataFlowGraphExecutionDrawer._determine_label(
                    computation, computation_label, status_abbreviation
                ),
                color=computation_color,
                tooltip=tooltip,
            )

        # Add data flows
        for (
            (from_computation, to_computation),
            data_flow_label,
        ) in data_flow_graph.edges_with_labels:
            inpud = computation_to_label_to_input[from_computation][data_flow_label]
            tooltip = "{label}: {inpud}".format(label=data_flow_label, inpud=inpud)
            # NOTE : the tooltip attribute for edge shows the defautl value '%3' at many edges alternative attributes
            # for tagging can be headtooltip, tailtooltip, as explained in http://www.graphviz.org/doc/info/attrs.html#d:tooltip
            digraph.edge(
                str(data_flow_graph.node_label(to_computation)),
                str(data_flow_graph.node_label(from_computation)),
                tooltip=tooltip,
            )

        digraph.render(format="svg")

    @staticmethod
    def _build_digraph() -> Digraph:
        return Digraph(
            directory=".",
            filename="data_flow_graph.gv",
            encoding="utf-8",
            engine="dot",
            node_attr={"color": "lightblue2", "style": "filled"},
            graph_attr={
                "packMode": "node",
                "margin": "0.2",
                "rankdir": "TB",
                "center": "true",
                "ratio": "0.7",
                "rank": "source",
                "nodesep": "0.1",
                "ranksep": "0.1",
                "splines": "true",
                "overlap": "vpsc",
            },
        )

    @staticmethod
    def _determine_status_abbreviation_and_color(
        computation: Computation,
        computation_to_process_with_status: Dict[
            Computation, Tuple[Process, _ProcessStatus]
        ]  # pylint: disable=invalid-name
    ) -> Tuple[str, str]:
        if isinstance(computation, File):
            return "F", "darkseagreen2"

        process, status = computation_to_process_with_status[computation]
        if status == _ProcessStatus.PENDING:
            return "P", "lightblue2"  # pending
        if status == _ProcessStatus.STARTED:
            if process.is_alive():
                return "R", "salmon2"  # running
            return "F", "darkseagreen2"  # finished
        raise NotImplementedError

    @staticmethod
    def _determine_label(
        computation: Computation,
        computation_label: ComputationLabel,
        status_abbreviation: str
    ) -> str:
        if isinstance(computation, File):
            return str(computation_label)

        return "{label}: {status}".format(
            label=computation_label, status=status_abbreviation
        )

    @staticmethod
    def _format_and_collect_inputs(
        computation: Computation, label_to_input: Dict[DataFlowLabel, Input]
    ) -> str:
        result = ""
        for input_label in computation.input_labels:
            inpud = computation.input(input_label)
            result += "%s : %s \n" % (input_label, str(inpud))
            label_to_input[input_label] = inpud
        return result

    @staticmethod
    def _format_output(output: Optional[Output]) -> str:
        if output is not None:
            return str(output) + "\n"
        return ""

    @staticmethod
    def _format_meta_data(data: Dict[str, Any]) -> str:
        result = ""
        for item_name, item_data in data.items():
            result += "%s : %s \n" % (item_name, str(item_data))
        return result
