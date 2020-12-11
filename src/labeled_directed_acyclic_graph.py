"""Node- and edge-labeled directed acyclic graphs."""

import copy
from abc import ABCMeta, abstractmethod
from typing import (
    Dict,
    Generator,
    Generic,
    ItemsView,
    Iterable,
    Iterator,
    KeysView,
    List,
    Set,
    Tuple,
    Type,
    TypeVar,
    ValuesView,
)

Node = TypeVar("Node")
"""Type variable for nodes."""  # pylint: disable=pointless-string-statement


NodeLabel = TypeVar("NodeLabel")
"""Type variable for node labels."""  # pylint: disable=pointless-string-statement


EdgeLabel = TypeVar("EdgeLabel")
"""Type variable for edge labels."""  # pylint: disable=pointless-string-statement


# As of python 3.8 we may use protocols instead.
# See https://www.python.org/dev/peps/pep-0544/
# and https://mypy.readthedocs.io/en/latest/protocols.html#simple-user-defined-protocols
class Equatable(metaclass=ABCMeta):
    """Interface for equatable nodes."""

    @abstractmethod
    def is_equivalent_to(self, other: object) -> bool:
        """Checks whether this node is equivalent to the other one.

        This method must be implemented by sub-classes in such a way that
        ``x.is_equivalent_to(y)`` is the same as ``y.is_equivalent_to(x)``, in
        other words, the relation ``is_equivalent_to`` is commutative.

        Args:
            other: The other node.
        """


"""Type variable for equatable nodes."""  # pylint: disable=pointless-string-statement
EquatableNode = TypeVar("EquatableNode", bound=Equatable)


# The graphs are acyclic by construction.
class LabeledDirectedAcyclicGraph(Generic[Node, NodeLabel, EdgeLabel]):
    """Representation of a generic node- and edge-labeled directed acyclic
    graph, commonly known by its acronym DAG, see
    https://en.wikipedia.org/wiki/Directed_acyclic_graph and
    https://en.wikipedia.org/wiki/Graph_labeling

    The graph is node- and edge-labeled, and acyclic by construction, in
    other words, you cannot create an instance of this class without those
    properties.

    This class is parameterized over the type variables ``Node``, ``NodeLabel``,
    and ``EdgeLabel``. Each instance of this class specifies a type for each
    of these type variables, for example, nodes could be ``int`` , node labels
    ``str`` and edge labels ``bool``, in which case a variable ``graph`` holding
    an instance of such a graph is declared by
    ``graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]``.

    Examples:
        Graphs can be constructed by multiple invocations of
        ``LabeledDirectedAcyclicGraph.add_node``:

        >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
        >>> graph.add_node("one", 1, {})
        >>> graph.add_node("two", 2, {})
        >>> graph.add_node("three", 3, {True: "one", False: "two"})
        >>> graph.add_node("four", 4, {True: "one", False: "three"})
        >>> graph.add_node("five", 5, {True: "one"})

        And they can be inspected using its various public properties and
        methods:

        >>> graph.roots
        [4, 5]
        >>> graph.nodes
        dict_keys([1, 2, 3, 4, 5])
        >>> graph.nodes_with_labels
        dict_items([(1, 'one'), (2, 'two'), (3, 'three'), (4, 'four'), (5, 'five')])
        >>> graph.edges
        dict_keys([(3, 1), (3, 2), (4, 1), (4, 3), (5, 1)])
        >>> graph.edges_with_labels
        dict_items([((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False), ((5, 1), True)])

        It is however hard to see the graph's structure in code from this way
        of construction. A more intuitive way is to construct graphs
        recursively by saying that the new graph consists of this and that
        graph joined by a new node, in other words, by merging existing graphs
        and connecting their roots with a new node. This is what the method
        ``LabeledDirectedAcyclicGraph.construct`` does. Let's see it in action.

        To merge graphs we need to be able to tell whether two nodes from
        different graphs are equivalent in which case they can be interchanged.
        To that end, a custom node type must inherit from the abstract class
        ``Equatable`` and implement the method ``Equatable.is_equivalent_to``, for
        example, as follows:

        >>> class Node(Equatable):
        ...     def __init__(self, value: str) -> None:
        ...         self.value = value
        ...
        ...     def is_equivalent_to(self, other: object) -> bool:
        ...         if type(other) is type(self):
        ...             return self.value == other.value
        ...         else:
        ...             return False
        ...
        ...     def __str__(self) -> str:
        ...         return "Node({value})".format(value=self.value)
        ...
        ...     def __repr__(self) -> str:
        ...         return str(self)

        As said, labeled, directed, and acyclic graphs with precisely one root
        can be constructed top-down with the method
        ``LabeledDirectedAcyclicGraph.construct``. It takes a node label, an
        equatable node that's going to be the root, and a map from edge labels
        to what's going to be subgraphs of the node.

        A singleton graph can be constructed as follows:

        >>> def singleton_graph() -> LabeledDirectedAcyclicGraph[Node, str, int]:
        ...     return LabeledDirectedAcyclicGraph.construct(
        ...         "singleton_root",
        ...         Node(value = "singleton_value"),
        ...         {},
        ...     )

        The graph's root is

        >>> graph = singleton_graph()
        >>> singleton_root = graph.root
        >>> str(singleton_root)
        'Node(singleton_value)'

        The root's label is

        >>> graph.node_label(singleton_root)
        'singleton_root'

        Based on existing graphs bigger and bigger graphs can be built by
        connecting the roots of the existing graphs with a new node that's
        going to be the root of the new graph. For example, based on the
        singleton graph a doubleton graph can be constructed as follows:

        >>> def doubleton_graph(suffix: str) -> LabeledDirectedAcyclicGraph[Node, str, int]:
        ...     return LabeledDirectedAcyclicGraph.construct(
        ...         "doubleton_root",
        ...         Node(value = "doubleton_value_{}".format(suffix)),
        ...         {1: singleton_graph()},
        ...     )

        The graph's root is labeled

        >>> doubleton_graph("x").root
        Node(doubleton_value_x)

        Not all labeled directed acyclic graphs can be constructed in this
        manner though, only exactly those that are weakly connected and have
        precisely one root.

        Let's extend our graph even more:

        >>> def tripleton_graph(doubleton_suffix: str, tripleton_suffix: str) -> LabeledDirectedAcyclicGraph[Node, str, int]:
        ...     return LabeledDirectedAcyclicGraph.construct(
        ...         "tripleton_root",
        ...         Node(value = "tripleton_value_{}".format(tripleton_suffix)),
        ...         {2: doubleton_graph(doubleton_suffix)},
        ...     )
        >>> def nton_graph(doubleton_suffix: str, tripleton_suffix: str) -> LabeledDirectedAcyclicGraph[Node, str, int]:
        ...     return LabeledDirectedAcyclicGraph.construct(
        ...         "nton_root",
        ...         Node(value = "nton_value"),
        ...         {
        ...             2: doubleton_graph(doubleton_suffix),
        ...             3: tripleton_graph(doubleton_suffix, tripleton_suffix),
        ...         },
        ...     )
        >>> str(nton_graph("x", "y").root)
        'Node(nton_value)'

        Note that all graphs given to ``LabeledDirectedAcyclicGraph.construct``
        must have at most one root. Under the hood it uses the static method
        ``LabeledDirectedAcyclicGraph.merge`` to merge the given graphs and the
        instance method ``LabeledDirectedAcyclicGraph.add_node`` to join the
        roots of the merged subgraphs by a new node.

        We wrap the construction of graphs inside functions to
        * make the construction parameterizable, for example, for data-flow
          graphs by different input-file paths and computation parameters;
        * make all parameters explicit in parameter lists, so that there is no
          need to search for local or global variables that are being used
          for parameterization somewhere in the code;
        * make the construction independent of the order of the definitions,
          for example, the function ``singleton_graph`` above could be defined
          after the function ``doubleton_graph`` without breaking the code
          (while this is not be advisable in the above case, there are graphs
          in which nodes do not have a linear order where all edges point in
          the same direction).

        A graph with multiple roots can be created as follows:

        >>> multiple_roots_graph = copy.copy(nton_graph("z", "z"))
        >>> multiple_roots_graph.add_node(
        ...     "value = another_root",
        ...     Node("another_value"),
        ...     {
        ...         2: "doubleton_root",
        ...         3: "tripleton_root",
        ...     },
        ... )
        >>> [str(root) for root in multiple_roots_graph.roots]
        ['Node(nton_value)', 'Node(another_value)']

        >>> multiple_roots_graph.root # doctest: +ELLIPSIS
        Traceback (most recent call last):
        AssertionError: The graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at ...>' has more than one root, namely 'Node(nton_value), Node(another_value)'

        A graph can be traversed depth first with the instance method
        ``LabeledDirectedAcyclicGraph.traverse_nodes_depth_first``.

        >>> graph = tripleton_graph("p", "q")
        >>> [str(node) for node in graph.traverse_nodes_depth_first(graph.root)]
        ['Node(singleton_value)', 'Node(doubleton_value_p)', 'Node(tripleton_value_q)']
    """

    @classmethod
    def format_edge(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        edge: Tuple[Node, Node],
    ) -> str:
        """Formats the given edge as string.

        Example:
            >>> LabeledDirectedAcyclicGraph.format_edge(('from', 'to'))
            '(from, to)'
        """
        return "({from_node}, {to_node})".format(from_node=edge[0], to_node=edge[1])

    @classmethod
    def construct(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        label: NodeLabel,
        node: EquatableNode,
        edge_label_to_graphs: Dict[
            EdgeLabel,
            "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
        ]
    ) -> "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]":
        """Constructs a new graph with the root ``node`` labeled ``label`` and
        edges labeled as given from the root to the roots of the given graphs.
        In other words, the given graphs are merged and their roots connected
        by the new node ``node``, which is then the root of the new graph.

        Note that each of the given graphs must have at most one root and
        that the returned graph has precisely one root.

        Args:
            label: The label of the root.
            node: The root.
            edge_label_to_graphs: A map from edge labels to graphs, where each
                edge label becomes the label of an edge from the root to the
                respective graph's root.

        Returns:
            A new graph with the new root that has the given sub-graphs.

        Raises:
            AssertionError: If at least one of the given graphs does not have
                at most one root.
            AssertionError: If the given graphs cannot be merged due to
                conflicts, namely,
                * if there is a node label that occurs in two graphs and the
                  corresponding nodes are not equivalent as defined by
                  ``Equatable.is_equivalent_to`` or have different numbers of
                  children;
                * if there is a node that has different labels or number of
                  children in different graphs;
                * if there is a node that has two edges with the same label in
                  two graphs to non-equivalent nodes as defined by
                  ``Equatable.is_equivalent_to``;
                * if there is an edge that has different labels in different
                  graphs.
                Note that instead of raising errors when two equivalent nodes
                have different numbers of children in different graphs, the
                children could be merged.

        Examples:
            For a step-by-step introduction on how to use this method see the
            class documentation of ``LabeledDirectedAcyclicGraph``.

            >>> class Node(Equatable):
            ...     def __init__(self, value: str) -> None:
            ...         self.value = value
            ...
            ...     def is_equivalent_to(self, other: object) -> bool:
            ...         if type(other) is type(self):
            ...             return self.value == other.value
            ...         else:
            ...             return False
            ...
            ...     def __str__(self) -> str:
            ...         return "Node({value})".format(value=self.value)
            ...
            ...     def __repr__(self) -> str:
            ...         return str(self)

            >>> singleton_graph = LabeledDirectedAcyclicGraph.construct(
            ...     "singleton_root",
            ...     Node(value = "singleton_value"),
            ...     {},
            ... )
            >>> (list(singleton_graph.nodes_with_labels), list(singleton_graph.edges_with_labels))
            ([(Node(singleton_value), 'singleton_root')], [])

            >>> doubleton_graph = LabeledDirectedAcyclicGraph.construct(
            ...     "doubleton_root",
            ...     Node(value = "doubleton_value"),
            ...     {"singleton": singleton_graph},
            ... )
            >>> list(doubleton_graph.nodes_with_labels)
            [(Node(singleton_value), 'singleton_root'), (Node(doubleton_value), 'doubleton_root')]
            >>> list(doubleton_graph.edges_with_labels)
            [((Node(doubleton_value), Node(singleton_value)), 'singleton')]

            >>> tripleton_graph = LabeledDirectedAcyclicGraph.construct(
            ...     "tripleton_root",
            ...     Node(value = "tripleton_value"),
            ...     {"doubleton": doubleton_graph},
            ... )
            >>> list(tripleton_graph.nodes_with_labels)
            [(Node(singleton_value), 'singleton_root'), (Node(doubleton_value), 'doubleton_root'), (Node(tripleton_value), 'tripleton_root')]
            >>> list(tripleton_graph.edges_with_labels)
            [((Node(doubleton_value), Node(singleton_value)), 'singleton'), ((Node(tripleton_value), Node(doubleton_value)), 'doubleton')]

            >>> nton_graph = LabeledDirectedAcyclicGraph.construct(
            ...     "nton_root",
            ...     Node(value = "nton_value"),
            ...     {
            ...         "doubleton": doubleton_graph,
            ...         "tripleton": tripleton_graph,
            ...     },
            ... )
            >>> list(nton_graph.nodes_with_labels)
            [(Node(singleton_value), 'singleton_root'), (Node(doubleton_value), 'doubleton_root'), (Node(tripleton_value), 'tripleton_root'), (Node(nton_value), 'nton_root')]
            >>> list(nton_graph.edges_with_labels)
            [((Node(doubleton_value), Node(singleton_value)), 'singleton'), ((Node(tripleton_value), Node(doubleton_value)), 'doubleton'), ((Node(nton_value), Node(doubleton_value)), 'doubleton'), ((Node(nton_value), Node(tripleton_value)), 'tripleton')]
        """

        result = LabeledDirectedAcyclicGraph.merge(edge_label_to_graphs.values())

        _edge_label_to_node = {}
        for _edge_label, _graph in edge_label_to_graphs.items():
            _edge_label_to_node[_edge_label] = _graph.node_label(_graph.root)

        result.add_node(label, node, _edge_label_to_node)
        return result

    @classmethod
    def merge(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        graphs: Iterable[
            "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
        ],
    ) -> "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]":
        """Merges the given graphs.

        Note that each of the given graphs must have at most one root.

        Args:
            graphs: The graphs to merge.

        Returns:
            The merge of the given graphs.

        Raises:
            AssertionError: If at least one of the given graphs does not have
                at most one root.
            AssertionError: If the given graphs cannot be merged due to
                conflicts, namely,
                * if there is a node label that occurs in two graphs and the
                  corresponding nodes are not equivalent as defined by
                  ``Equatable.is_equivalent_to`` or have different numbers of
                  children;
                * if there is a node that has different labels or number of
                  children in different graphs;
                * if there is a node that has two edges with the same label in
                  two graphs to non-equivalent nodes as defined by
                  ``Equatable.is_equivalent_to``;
                * if there is an edge that has different labels in different
                  graphs.
                Note that instead of raising errors when two equivalent nodes
                have different numbers of children in different graphs, the
                children could be merged.

        Examples:
            First, define a custom equatable-node type.

            >>> class Node(Equatable):
            ...     def __init__(self, value: int) -> None:
            ...         self.value = value
            ...
            ...     def is_equivalent_to(self, other: object) -> bool:
            ...         if type(other) is type(self):
            ...             return self.value == other.value
            ...         else:
            ...             return False
            ...
            ...     def __str__(self) -> str:
            ...         return "Node({value})".format(value=self.value)
            ...
            ...     def __repr__(self) -> str:
            ...         return str(self)

            Merge no graphs

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([], [])

            Merge one empty graph

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([LabeledDirectedAcyclicGraph()])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([], [])

            Merge three empty graphs

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([LabeledDirectedAcyclicGraph(), LabeledDirectedAcyclicGraph(), LabeledDirectedAcyclicGraph()])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([], [])

            Merge one singleton graph

            >>> singleton_graph = LabeledDirectedAcyclicGraph()
            >>> singleton_graph.add_node("a", Node(value = 1), {})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([singleton_graph])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([(Node(1), 'a')], [])

            Merge one singleton graph and two empty graphs

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([singleton_graph, LabeledDirectedAcyclicGraph(), LabeledDirectedAcyclicGraph()])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([(Node(1), 'a')], [])

            Merge one singleton graph, a structurally identical one and an empty graph

            >>> copied_singleton_graph = copy.copy(singleton_graph)
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([singleton_graph, copied_singleton_graph, LabeledDirectedAcyclicGraph()])
            >>> (list(merged_graph.nodes_with_labels), list(merged_graph.edges_with_labels))
            ([(Node(1), 'a')], [])

            Merge one ``n``-ton graph

            >>> nton_graph = LabeledDirectedAcyclicGraph()
            >>> nton_graph.add_node("a", Node(value = 1), {})
            >>> nton_graph.add_node("b", Node(value = 2), {})
            >>> nton_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> nton_graph.add_node("d", Node(value = 4), {True: "a", False: "c"})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False)]

            Merge one ``n``-ton graph and one empty graph

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, LabeledDirectedAcyclicGraph()])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False)]

            Merge one ``n``-ton graph and one singleton graph with an overlapping node

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, singleton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False)]

            Merge one ``n``-ton graph and one singleton graph with another node

            >>> another_singleton_graph = LabeledDirectedAcyclicGraph()
            >>> another_singleton_graph.add_node("q", Node(value = 0), {})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, another_singleton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd'), (Node(0), 'q')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False)]

            Merge one ``n``-ton graph, a structurally identical one and an empty graph

            >>> copied_nton_graph = copy.copy(nton_graph)
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, copied_nton_graph, LabeledDirectedAcyclicGraph()])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False)]

            Merge two non-overlapping ``n``-ton graphs

            >>> non_overlapping_nton_graph = LabeledDirectedAcyclicGraph()
            >>> non_overlapping_nton_graph.add_node("w", Node(value = 11), {})
            >>> non_overlapping_nton_graph.add_node("x", Node(value = 12), {False: "w"})
            >>> non_overlapping_nton_graph.add_node("y", Node(value = 13), {False: "x"})
            >>> non_overlapping_nton_graph.add_node("z", Node(value = 14), {True: "w", False: "y"})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, non_overlapping_nton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd'), (Node(11), 'w'), (Node(12), 'x'), (Node(13), 'y'), (Node(14), 'z')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False), ((Node(12), Node(11)), False), ((Node(13), Node(12)), False), ((Node(14), Node(11)), True), ((Node(14), Node(13)), False)]

            Merge two partially overlapping ``n``-ton graphs

            >>> partially_overlapping_nton_graph = LabeledDirectedAcyclicGraph()
            >>> partially_overlapping_nton_graph.add_node("b", Node(value = 2), {})
            >>> partially_overlapping_nton_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> partially_overlapping_nton_graph.add_node("x", Node(value = 12), {False: "c"})
            >>> partially_overlapping_nton_graph.add_node("y", Node(value = 13), {False: "x"})
            >>> partially_overlapping_nton_graph.add_node("z", Node(value = 14), {True: "x", False: "y"})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, partially_overlapping_nton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd'), (Node(12), 'x'), (Node(13), 'y'), (Node(14), 'z')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False), ((Node(12), Node(3)), False), ((Node(13), Node(12)), False), ((Node(14), Node(12)), True), ((Node(14), Node(13)), False)]

            Merge one ``n``-ton graph and one that includes the other

            >>> subsuming_nton_graph = LabeledDirectedAcyclicGraph()
            >>> subsuming_nton_graph.add_node("a", Node(value = 1), {})
            >>> subsuming_nton_graph.add_node("b", Node(value = 2), {})
            >>> subsuming_nton_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> subsuming_nton_graph.add_node("d", Node(value = 4), {True: "a", False: "c"})
            >>> subsuming_nton_graph.add_node("x", Node(value = 12), {False: "c"})
            >>> subsuming_nton_graph.add_node("y", Node(value = 13), {False: "x", True: "d"})
            >>> subsuming_nton_graph.add_node("z", Node(value = 14), {True: "x", False: "y"})
            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, subsuming_nton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd'), (Node(12), 'x'), (Node(13), 'y'), (Node(14), 'z')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False), ((Node(12), Node(3)), False), ((Node(13), Node(12)), False), ((Node(13), Node(4)), True), ((Node(14), Node(12)), True), ((Node(14), Node(13)), False)]

            Merge multiple ``n``-ton and singleton graphs

            >>> merged_graph = LabeledDirectedAcyclicGraph.merge([nton_graph, copied_nton_graph, non_overlapping_nton_graph, singleton_graph, copied_singleton_graph])
            >>> list(merged_graph.nodes_with_labels)
            [(Node(1), 'a'), (Node(2), 'b'), (Node(3), 'c'), (Node(4), 'd'), (Node(11), 'w'), (Node(12), 'x'), (Node(13), 'y'), (Node(14), 'z')]
            >>> list(merged_graph.edges_with_labels)
            [((Node(3), Node(2)), False), ((Node(4), Node(1)), True), ((Node(4), Node(3)), False), ((Node(12), Node(11)), False), ((Node(13), Node(12)), False), ((Node(14), Node(11)), True), ((Node(14), Node(13)), False)]

            Try to merge graphs where at least one does not have at most one root

            >>> multiple_roots_graph = LabeledDirectedAcyclicGraph()
            >>> multiple_roots_graph.add_node("a", Node(value = 1), {})
            >>> multiple_roots_graph.add_node("b", Node(value = 2), {})
            >>> multiple_roots_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, multiple_roots_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' has more than one root, namely 'Node(1), Node(3)'

            Try to merge graphs with conflicts.
            First, there is a node label that occurs in two graphs and the
            corresponding nodes are not equivalent:

            >>> conflicting_node_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_node_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_node_graph.add_node("b", Node(value = 2), {})
            >>> conflicting_node_graph.add_node("c", Node(value = 42), {False: "b"})
            >>> conflicting_node_graph.add_node("d", Node(value = 4), {True: "a", False: "c"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_node_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The label 'c' has non-equivalent nodes in the two graphs, namely 'Node(3)' and 'Node(42)'

            Secondly, two nodes with the same label have different numbers of
            children:

            >>> conflicting_node_label_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_node_label_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_node_label_graph.add_node("b", Node(value = 2), {})
            >>> conflicting_node_label_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> conflicting_node_label_graph.add_node("d", nton_graph.node("d"), {True: "a"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_node_label_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The same node 'Node(4)' with the labels 'd' and 'd' has different numbers of children in the two graphs, namely '2' and '1'
            * The nodes 'Node(4)' and 'Node(4)' with the same label 'd' in the two graphs have different numbers of children, namely '2' and '1'

            Thirdly, there is a node that has different labels in different
            graphs:

            >>> conflicting_node_label_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_node_label_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_node_label_graph.add_node("b", Node(value = 2), {})
            >>> conflicting_node_label_graph.add_node("q", nton_graph.node("c"), {False: "b"})
            >>> conflicting_node_label_graph.add_node("d", Node(value = 4), {True: "a", False: "q"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_node_label_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The node 'Node(3)' has different labels in the two graphs, namely 'c' and 'q'

            Fourthly, two nodes with the same label have different numbers of
            children:

            >>> conflicting_node_label_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_node_label_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_node_label_graph.add_node("b", Node(value = 2), {})
            >>> conflicting_node_label_graph.add_node("c", Node(value = 3), {False: "b"})
            >>> conflicting_node_label_graph.add_node("d", Node(value = 4), {True: "a"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_node_label_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The nodes 'Node(4)' and 'Node(4)' with the same label 'd' in the two graphs have different numbers of children, namely '2' and '1'

            Fifthly, there is a node that has two edges with the same label in
            two graphs to non-equivalent nodes:

            >>> conflicting_edge_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_edge_graph.add_node("z", Node(value = 0), {})
            >>> conflicting_edge_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_edge_graph.add_node("c", Node(value = 3), {False: "z"})
            >>> conflicting_edge_graph.add_node("d", Node(value = 4), {True: "a", False: "c"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_edge_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The edge from node 'Node(3)' with label 'False' has different target nodes in the two graphs, namely 'Node(2)' and 'Node(0)'

            Sixthly, there is an edge that has different labels in
            different graphs:

            >>> conflicting_edge_label_graph = LabeledDirectedAcyclicGraph()
            >>> conflicting_edge_label_graph.add_node("a", Node(value = 1), {})
            >>> conflicting_edge_label_graph.add_node("b", Node(value = 2), {})
            >>> conflicting_edge_label_graph.add_node("c", Node(value = 3), {True: "b"})
            >>> conflicting_edge_label_graph.add_node("d", Node(value = 4), {True: "a", False: "c"})
            >>> LabeledDirectedAcyclicGraph.merge([nton_graph, conflicting_edge_label_graph, partially_overlapping_nton_graph]) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The two graphs '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' and '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:
            * The edge '(Node(3), Node(2))' has different labels in the two graphs, namely 'False' and 'True'
        """

        result = (
            LabeledDirectedAcyclicGraph()
        )  # type: LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]
        for _graph in graphs:
            if _graph:
                merge_conflicts = LabeledDirectedAcyclicGraph._merge_conflicts(
                    result, _graph
                )
                assert not merge_conflicts, (
                    "The two graphs '"
                    + str(result)
                    + "' and '"
                    + str(_graph)
                    + "' cannot be merged due to conflicting nodes and/or edges and their labels, namely due to the following conflicts:\n* "
                    + "\n* ".join(merge_conflicts)
                )
                for _graph_node in _graph.traverse_nodes_depth_first(_graph.root):
                    if not result.nodes:
                        result.add_node(
                            _graph.node_label(_graph_node),
                            _graph_node,
                            _graph.edge_label_to_child_node_label(_graph_node),
                        )
                    else:
                        if not (
                            all(
                                _graph_node.is_equivalent_to(result_node)
                                for result_node in result.nodes
                            )
                            or (result.has_node_label(_graph.node_label(_graph_node)))
                        ):
                            result.add_node(
                                _graph.node_label(_graph_node),
                                _graph_node,
                                _graph.edge_label_to_child_node_label(_graph_node),
                            )
        return result

    @classmethod
    def _merge_conflicts(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        first_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
        second_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
    ) -> List[str]:
        nodes_conflicts = LabeledDirectedAcyclicGraph._nodes_merge_conflicts(
            first_graph, second_graph
        )
        if nodes_conflicts:
            return nodes_conflicts
        return LabeledDirectedAcyclicGraph._edges_merge_conflicts(
            first_graph, second_graph
        )

    @classmethod
    def _nodes_merge_conflicts(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        first_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
        second_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"
    ) -> List[str]:
        conflicts = []
        for first_node, first_label in first_graph.nodes_with_labels:
            if second_graph.has_node(first_node):
                second_label = second_graph.node_label(first_node)
                if second_label != first_label:
                    conflicts.append(
                        "The node '"
                        + str(first_node)
                        + "' has different labels in the two graphs, namely '"
                        + str(first_label)
                        + "' and '"
                        + str(second_label)
                        + "'"
                    )
                first_children_count = len(first_graph.child_nodes(first_node))
                second_children_count = len(second_graph.child_nodes(first_node))
                if first_children_count != second_children_count:
                    conflicts.append(
                        "The same node '"
                        + str(first_node)
                        + "' with the labels '"
                        + str(first_label)
                        + "' and '"
                        + str(second_label)
                        + "' has different numbers of children in the two graphs, namely '"
                        + str(first_children_count)
                        + "' and '"
                        + str(second_children_count)
                        + "'"
                    )
            if second_graph.has_node_label(first_label):
                second_node = second_graph.node(first_label)
                if not second_node.is_equivalent_to(first_node):
                    conflicts.append(
                        "The label '"
                        + str(first_label)
                        + "' has non-equivalent nodes in the two graphs, namely '"
                        + str(first_node)
                        + "' and '"
                        + str(second_node)
                        + "'"
                    )
                first_children_count = len(first_graph.child_nodes(first_node))
                second_children_count = len(second_graph.child_nodes(second_node))
                if first_children_count != second_children_count:
                    conflicts.append(
                        "The nodes '"
                        + str(first_node)
                        + "' and '"
                        + str(second_node)
                        + "' with the same label '"
                        + str(first_label)
                        + "' in the two graphs have different numbers of children, namely '"
                        + str(first_children_count)
                        + "' and '"
                        + str(second_children_count)
                        + "'"
                    )
        return conflicts

    @classmethod
    def _edges_merge_conflicts(
        cls: Type["LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"],
        first_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]",
        second_graph: "LabeledDirectedAcyclicGraph[EquatableNode, NodeLabel, EdgeLabel]"
    ) -> List[str]:
        conflicts = []
        for first_edge, first_edge_label in first_graph.edges_with_labels:
            first_from_node, first_to_node = first_edge
            from_node_label = first_graph.node_label(first_from_node)
            to_node_label = first_graph.node_label(first_to_node)
            if second_graph.has_node_label(from_node_label):
                second_from_node = second_graph.node(from_node_label)
                if second_graph.has_node_label(to_node_label):
                    second_to_node = second_graph.node(to_node_label)
                    second_edge = (second_from_node, second_to_node)
                    if second_graph.has_edge(second_edge):
                        second_edge_label = second_graph.edge_label(second_edge)
                        if second_edge_label != first_edge_label:
                            conflicts.append(
                                "The edge '"
                                + LabeledDirectedAcyclicGraph.format_edge(first_edge)
                                + "' has different labels in the two graphs, namely '"
                                + str(first_edge_label)
                                + "' and '"
                                + str(second_edge_label)
                                + "'"
                            )
                if second_graph.has_edge_label(second_from_node, first_edge_label):
                    second_to_node = second_graph.child_node(
                        second_from_node, first_edge_label
                    )
                    if not first_to_node.is_equivalent_to(second_to_node):
                        conflicts.append(
                            "The edge from node '"
                            + str(first_from_node)
                            + "' with label '"
                            + str(first_edge_label)
                            + "' has different target nodes in the two graphs, namely '"
                            + str(first_to_node)
                            + "' and '"
                            + str(second_to_node)
                            + "'"
                        )
        return conflicts

    def __init__(self) -> None:
        """The constructor initializes an empty graph.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph()
            >>> (list(graph.nodes_with_labels), list(graph.edges_with_labels))
            ([], [])
            >>> len(graph)
            0
            >>> 42 in graph
            False
            >>> [node for node in graph]
            []

            >>> copied_graph = copy.copy(graph)
            >>> (list(copied_graph.nodes_with_labels), list(copied_graph.edges_with_labels))
            ([], [])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> list(graph.nodes_with_labels)
            [(0, 'zero'), (1, 'one'), (2, 'two'), (3, 'three'), (4, 'four'), (5, 'five')]
            >>> list(graph.edges_with_labels)
            [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False), ((5, 1), True)]

            >>> len(graph)
            6
            >>> 0 in graph and 1 in graph and 2 in graph and 3 in graph and 4 in graph and 5 in graph
            True
            >>> 42 in graph
            False
            >>> [node for node in graph]
            [0, 1, 2, 3, 4, 5]

            >>> copied_graph = copy.copy(graph)
            >>> list(copied_graph.nodes_with_labels)
            [(0, 'zero'), (1, 'one'), (2, 'two'), (3, 'three'), (4, 'four'), (5, 'five')]
            >>> list(copied_graph.edges_with_labels)
            [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False), ((5, 1), True)]
        """
        self._node_to_label = {}  # type: Dict[Node, NodeLabel]
        self._label_to_node = {}  # type: Dict[NodeLabel, Node]
        self._edge_to_label = {}  # type: Dict[Tuple[Node, Node], EdgeLabel]
        self._node_to_child_nodes = {}  # type: Dict[Node, List[Node]]
        self._node_to_parent_nodes = {}  # type: Dict[Node, List[Node]]
        self._node_to_label_to_child_node = (
            {}
        )  # type: Dict[Node, Dict[EdgeLabel, Node]]

    def __len__(self) -> int:
        """Makes ``len(graph)`` return the number of nodes."""
        return len(self._node_to_label)

    def __contains__(self, node: Node) -> bool:
        """Makes ``node in graph`` tell whether ``node`` is contained in ``graph``."""
        return node in self._node_to_label

    def __iter__(self) -> Iterator[Node]:
        """Makes ``for node in graph`` iterate over all nodes in ``graph``."""
        return self.nodes.__iter__()

    def __copy__(self) -> "LabeledDirectedAcyclicGraph[Node, NodeLabel, EdgeLabel]":
        """Makes ``copy.copy(graph)`` make a shallow copy of ``graph``."""
        copied_graph = (
            LabeledDirectedAcyclicGraph()
        )  # type: LabeledDirectedAcyclicGraph[Node, NodeLabel, EdgeLabel]
        copied_graph._node_to_label = copy.copy(self._node_to_label)
        copied_graph._label_to_node = copy.copy(self._label_to_node)
        copied_graph._edge_to_label = copy.copy(self._edge_to_label)
        copied_graph._node_to_child_nodes = copy.copy(self._node_to_child_nodes)
        copied_graph._node_to_parent_nodes = copy.copy(self._node_to_parent_nodes)
        copied_graph._node_to_label_to_child_node = copy.copy(
            self._node_to_label_to_child_node
        )
        return copied_graph

    @property
    def nodes(self) -> KeysView[Node]:
        """All nodes

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.nodes
            dict_keys([])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.nodes
            dict_keys([0, 1, 2, 3, 4, 5])
        """
        return self._node_to_label.keys()

    @property
    def nodes_with_labels(self) -> ItemsView[Node, NodeLabel]:
        """All pairs of nodes and labels

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.nodes_with_labels
            dict_items([])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.nodes_with_labels
            dict_items([(0, 'zero'), (1, 'one'), (2, 'two'), (3, 'three'), (4, 'four'), (5, 'five')])
        """
        return self._node_to_label.items()

    @property
    def node_labels(self) -> ValuesView[NodeLabel]:
        """All node labels

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.node_labels
            dict_values([])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.node_labels
            dict_values(['zero', 'one', 'two', 'three', 'four', 'five'])
        """
        return self._node_to_label.values()

    @property
    def edges(self) -> KeysView[Tuple[Node, Node]]:
        """All edges

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.edges
            dict_keys([(3, 1), (3, 2), (4, 1), (4, 3)])
        """
        return self._edge_to_label.keys()

    @property
    def edges_with_labels(self) -> ItemsView[Tuple[Node, Node], EdgeLabel]:
        """All pairs of edges and labels

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.edges_with_labels
            dict_items([])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.edges_with_labels
            dict_items([((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False), ((5, 1), True)])
        """
        return self._edge_to_label.items()

    @property
    def edge_labels(self) -> ValuesView[EdgeLabel]:
        """All edge labels

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.edge_labels
            dict_values([])

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.edge_labels
            dict_values([True, False, True, False, True])
        """
        return self._edge_to_label.values()

    @property
    def roots(self) -> List[Node]:
        """All roots

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.roots
            []

            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.roots
            [4]

            >>> graph.add_node("zero", 0, {})
            >>> graph.roots
            [4, 0]

            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.roots
            [4, 0, 5]
        """
        return [node for node in self.nodes if not self.parent_nodes(node)]

    @property
    def root(self) -> Node:
        """The root

        Raises:
            AssertionError: If there is no root or there is more than 1 root.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.root # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' does not have a root

            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})

            >>> graph.root
            4

            >>> graph.add_node("zero", 0, {})
            >>> graph.root # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' has more than one root, namely '4, 0'

            >>> graph.add_node("five", 5, {True: "one"})
            >>> graph.root # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>' has more than one root, namely '4, 0, 5'
        """
        roots = self.roots
        assert len(roots) != 0, "The graph '" + str(self) + "' does not have a root"
        assert len(roots) == 1, (
            "The graph '"
            + str(self)
            + "' has more than one root, namely '"
            + ", ".join([str(root) for root in roots])
            + "'"
        )
        return roots[0]

    def add_node(
        self,
        label: NodeLabel,
        node: Node,
        edge_label_to_child_node_label: Dict[EdgeLabel, NodeLabel],
    ) -> None:
        """Adds the node ``node`` with label ``label`` and labeled edges
        to child nodes as specified by ``edge_label_to_child_node_label``.

        Note that nodes for the given child node labels must exist.

        Args:
            label: The node's label.
            node: The node itself.
            edge_label_to_child_node_label: The node's outgoing labeled edges
                to child nodes, where the latter are identified by their
                labels.

        Raises:
            AssertionError: If the node ``node`` already exists.
            AssertionError: If the label ``label`` already exists.
            AssertionError: If one of the given child node labels does not
                exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> list(graph.nodes_with_labels)
            [(0, 'zero'), (1, 'one'), (2, 'two'), (3, 'three'), (4, 'four'), (5, 'five')]

            >>> list(graph.edges_with_labels)
            [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False), ((5, 1), True), ((5, 2), False)]

            >>> graph.add_node("four", 42, {}) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: A node, namely '4', with label 'four' does already exist in the graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'

            >>> graph.add_node("seven", 4, {}) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The node '4' already exists in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'

            >>> graph.add_node("seven", 7, {True: "fourty-two"}) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: There is no node with label 'fourty-two' in the graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        assert node not in self._node_to_label, (
            "The node '"
            + str(node)
            + "' already exists in the directed acyclic graph '"
            + str(self)
            + "'"
        )
        assert label not in self._label_to_node, (
            "A node, namely '"
            + str(self._label_to_node[label])
            + "', with label '"
            + str(label)
            + "' does already exist in the graph '"
            + str(self)
            + "'"
        )
        self._node_to_label[node] = label
        self._label_to_node[label] = node
        self._node_to_child_nodes[node] = []
        self._node_to_parent_nodes[node] = []
        if node not in self._node_to_label_to_child_node:
            self._node_to_label_to_child_node[node] = {}

        for edge_label, child_node_label in edge_label_to_child_node_label.items():
            child_node = self.node(child_node_label)
            self._node_to_child_nodes[node].append(child_node)
            self._node_to_parent_nodes[child_node].append(node)
            self._edge_to_label[(node, child_node)] = edge_label
            self._node_to_label_to_child_node[node][edge_label] = child_node

    def has_node(self, node: Node) -> bool:
        """Checks whether the node ``node`` exists.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.has_node(42)
            False

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.has_node(0)
            True

            >>> graph.has_node(1)
            True

            >>> graph.has_node(2)
            True

            >>> graph.has_node(3)
            True

            >>> graph.has_node(4)
            True

            >>> graph.has_node(5)
            True

            >>> graph.has_node(42)
            False
        """
        return node in self._node_to_label

    def node(self, label: NodeLabel) -> Node:
        """Returns the node with the label ``label``.

        Raises:
            AssertionError: If the label ``label`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.node("zero")
            0

            >>> graph.node("one")
            1

            >>> graph.node("two")
            2

            >>> graph.node("three")
            3

            >>> graph.node("four")
            4

            >>> graph.node("five")
            5

            >>> graph.node("fourty-two") # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: There is no node with label 'fourty-two' in the graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node_label(label)
        return self._label_to_node[label]

    def has_node_label(self, label: NodeLabel) -> bool:
        """Checks whether the node label ``label`` exists.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.has_node_label("fourty-two")
            False

            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.has_node_label("zero")
            True

            >>> graph.has_node_label("one")
            True

            >>> graph.has_node_label("two")
            True

            >>> graph.has_node_label("three")
            True

            >>> graph.has_node_label("four")
            True

            >>> graph.has_node_label("five")
            True

            >>> graph.has_node_label("fourty-two")
            False
        """
        return label in self._label_to_node

    def node_label(self, node: Node) -> NodeLabel:
        """Returns the label of the node ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.node_label(0)
            'zero'

            >>> graph.node_label(1)
            'one'

            >>> graph.node_label(2)
            'two'

            >>> graph.node_label(3)
            'three'

            >>> graph.node_label(4)
            'four'

            >>> graph.node_label(5)
            'five'

            >>> graph.node_label(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return self._node_to_label[node]

    def child_nodes(self, node: Node) -> List[Node]:
        """Returns the child nodes of the node ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})

            >>> graph.child_nodes(0)
            []

            >>> graph.child_nodes(1)
            []

            >>> graph.child_nodes(2)
            []

            >>> graph.child_nodes(3)
            [1, 2]

            >>> graph.child_nodes(4)
            [1, 3]

            >>> graph.child_nodes(5)
            [1]

            >>> graph.child_nodes(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return self._node_to_child_nodes[node]

    def child_nodes_with_labels(
        self, node: Node
    ) -> List[Tuple[Node, NodeLabel, EdgeLabel]]:
        """Returns a list of triples with one element for each child node of the
        node ``node``, where the first component is the child node itself, the
        second its label, and the third the label of the edge from the node
        ``node`` to the child node.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})

            >>> graph.child_nodes_with_labels(0)
            []

            >>> graph.child_nodes_with_labels(1)
            []

            >>> graph.child_nodes_with_labels(2)
            []

            >>> graph.child_nodes_with_labels(3)
            [(1, 'one', True), (2, 'two', False)]

            >>> graph.child_nodes_with_labels(4)
            [(1, 'one', True), (3, 'three', False)]

            >>> graph.child_nodes_with_labels(5)
            [(1, 'one', True)]

            >>> graph.child_nodes_with_labels(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return [
            (
                child_node,
                self._node_to_label[child_node],
                self._edge_to_label[(node, child_node)],
            )
            for child_node in self._node_to_child_nodes[node]
        ]

    def has_child_node(self, node: Node, label: EdgeLabel) -> bool:
        """Checks whether the node ``node`` has an edge with label ``label`` to
        a child node.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})

            >>> graph.has_child_node(0, True)
            False

            >>> graph.has_child_node(0, False)
            False

            >>> graph.has_child_node(1, True)
            False

            >>> graph.has_child_node(1, False)
            False

            >>> graph.has_child_node(3, True)
            True

            >>> graph.has_child_node(3, False)
            True

            >>> graph.has_child_node(5, True)
            True

            >>> graph.has_child_node(5, False)
            False

            >>> graph.has_child_node(42, False) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return label in self._node_to_label_to_child_node[node]

    def child_node(self, node: Node, label: EdgeLabel) -> Node:
        """Returns the child node of the node ``node`` that is the target of the
        edge with label ``label`` from the source ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.
            AssertionError: If there is no edge with label ``label`` from the
                source ``node``.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one"})

            >>> graph.child_node(3, True)
            1

            >>> graph.child_node(3, False)
            2

            >>> graph.child_node(4, True)
            1

            >>> graph.child_node(4, False)
            3

            >>> graph.child_node(5, True)
            1

            >>> graph.child_node(42, False) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'

            >>> graph.child_node(1, True) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: There is no edge from the node '1' with label 'True' in the graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        self._assert_existence_of_edge_label(node, label)
        return self._node_to_label_to_child_node[node][label]

    def parent_nodes(self, node: Node) -> List[Node]:
        """Returns the parent nodes of the node ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.parent_nodes(0)
            []

            >>> graph.parent_nodes(1)
            [3, 4, 5]

            >>> graph.parent_nodes(2)
            [3, 5]

            >>> graph.parent_nodes(3)
            [4]

            >>> graph.parent_nodes(4)
            []

            >>> graph.parent_nodes(5)
            []

            >>> graph.parent_nodes(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return self._node_to_parent_nodes[node]

    def parent_nodes_with_labels(
        self, node: Node
    ) -> List[Tuple[Node, NodeLabel, EdgeLabel]]:
        """Returns a list of triples with one element for each parent node of the
        node ``node``, where the first component is the parent node itself, the
        second its label, and the third the label of the edge from the parent node
        to the node ``node``.

        Note that the edge labels of the parent nodes of the node ``node`` are
        in general not unique.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.parent_nodes_with_labels(0)
            []

            >>> graph.parent_nodes_with_labels(1)
            [(3, 'three', True), (4, 'four', True), (5, 'five', True)]

            >>> graph.parent_nodes_with_labels(2)
            [(3, 'three', False), (5, 'five', False)]

            >>> graph.parent_nodes_with_labels(3)
            [(4, 'four', False)]

            >>> graph.parent_nodes_with_labels(4)
            []

            >>> graph.parent_nodes_with_labels(5)
            []

            >>> graph.parent_nodes_with_labels(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return [
            (
                parent_node,
                self._node_to_label[parent_node],
                self._edge_to_label[(parent_node, node)],
            )
            for parent_node in self._node_to_parent_nodes[node]
        ]

    def ancestors(self, node: Node) -> Set[Node]:
        """Returns the ancestors of the node ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.ancestors(0)
            set()

            >>> graph.ancestors(1)
            {3, 4, 5}

            >>> graph.ancestors(2)
            {3, 4, 5}

            >>> graph.ancestors(3)
            {4}

            >>> graph.ancestors(4)
            set()

            >>> graph.ancestors(5)
            set()

            >>> graph.ancestors(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        result = set()
        nodes_to_visit = [node]
        while nodes_to_visit:
            current_node = nodes_to_visit.pop()
            for parent_node in self.parent_nodes(current_node):
                nodes_to_visit.append(parent_node)
                result.add(parent_node)
        return result

    def has_edge(self, edge: Tuple[Node, Node]) -> bool:
        """Checks whether the edge ``edge`` exists.

        Raises:
            AssertionError: If the source or target node of the edge ``edge``
                does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.has_edge((1, 1)) or graph.has_edge((1, 2)) or graph.has_edge((1, 3)) or graph.has_edge((1, 4))
            False

            >>> graph.has_edge((2, 1)) or graph.has_edge((2, 2)) or graph.has_edge((2, 3)) or graph.has_edge((2, 4))
            False

            >>> graph.has_edge((3, 1)) and graph.has_edge((3, 2)) and not graph.has_edge((3, 3)) and not graph.has_edge((3, 4))
            True

            >>> graph.has_edge((4, 1)) and not graph.has_edge((4, 2)) and graph.has_edge((4, 3)) and not graph.has_edge((4, 4))
            True

            >>> graph.has_edge((1, 42)) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'

            >>> graph.has_edge((42, 1)) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(edge[0])
        self._assert_existence_of_node(edge[1])
        return edge in self._edge_to_label

    def edge(self, node: Node, label: EdgeLabel) -> Tuple[Node, Node]:
        """Returns the edge from the node ``node`` with the label ``label``.

        Raises:
            AssertionError: If the node ``node`` does not exist.
            AssertionError: If there is no edge from the node ``node`` with the
                label ``label``.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.edge(3, True)
            (3, 1)

            >>> graph.edge(3, False)
            (3, 2)

            >>> graph.edge(4, True)
            (4, 1)

            >>> graph.edge(4, False)
            (4, 3)

            >>> graph.edge(42, True) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'

            >>> graph.edge(4, None) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: There is no edge from the node '4' with label 'None' in the graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        return (node, self.child_node(node, label))

    def has_edge_label(self, node: Node, label: EdgeLabel) -> bool:
        """Checks whether the node ``node`` has an edge with label ``label``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.has_edge_label(1, True)
            False

            >>> graph.has_edge_label(1, False)
            False

            >>> graph.has_edge_label(3, True)
            True

            >>> graph.has_edge_label(3, False)
            True

            >>> graph.has_edge_label(42, True) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return label in self._node_to_label_to_child_node[node]

    def edge_label(self, edge: Tuple[Node, Node]) -> EdgeLabel:
        """Returns the label of the edge ``edge``.

        Raises:
            AssertionError: If the edge ``edge`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.edge_label((3, 1))
            True

            >>> graph.edge_label((3, 2))
            False

            >>> graph.edge_label((4, 2)) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '(4, 2)' is not an edge in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_edge(edge)
        return self._edge_to_label[edge]

    def edge_label_to_child_node_label(self, node: Node) -> Dict[EdgeLabel, NodeLabel]:
        """Returns the map from edge labels to child node labels for edges (and
        child nodes) of the node ``node``

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> graph.edge_label_to_child_node_label(0)
            {}

            >>> graph.edge_label_to_child_node_label(1)
            {}

            >>> graph.edge_label_to_child_node_label(2)
            {}

            >>> graph.edge_label_to_child_node_label(3)
            {True: 'one', False: 'two'}

            >>> graph.edge_label_to_child_node_label(4)
            {True: 'one', False: 'three'}

            >>> graph.edge_label_to_child_node_label(5)
            {True: 'one', False: 'two'}

            >>> graph.edge_label_to_child_node_label(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        return {
            self.edge_label((node, child_node)): self.node_label(child_node)
            for child_node in self.child_nodes(node)
        }

    def traverse_leaves(self) -> Generator[Node, None, None]:
        """Traverses the leaf nodes.

        Yields:
            One leaf after another.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> list(graph.traverse_leaves())
            []

            >>> graph.add_node("zero", 0, {})
            >>> list(graph.traverse_leaves())
            [0]

            >>> graph.add_node("one", 1, {})
            >>> list(graph.traverse_leaves())
            [0, 1]

            >>> graph.add_node("two", 2, {})
            >>> list(graph.traverse_leaves())
            [0, 1, 2]

            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> list(graph.traverse_leaves())
            [0, 1, 2]

            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> list(graph.traverse_leaves())
            [0, 1, 2]

            >>> graph.add_node("five", 5, {True: "one", False: "two"})
            >>> list(graph.traverse_leaves())
            [0, 1, 2]
        """
        for node in self.nodes:
            if not self.child_nodes(node):
                yield node

    def traverse_roots(self) -> Generator[Node, None, None]:
        """Traverses the root nodes.

        Yields:
            One root after another.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> list(graph.traverse_roots())
            []

            >>> graph.add_node("zero", 0, {})
            >>> list(graph.traverse_roots())
            [0]

            >>> graph.add_node("one", 1, {})
            >>> list(graph.traverse_roots())
            [0, 1]

            >>> graph.add_node("two", 2, {})
            >>> list(graph.traverse_roots())
            [0, 1, 2]

            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> list(graph.traverse_roots())
            [0, 3]

            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> list(graph.traverse_roots())
            [0, 4]

            >>> graph.add_node("five", 5, {True: "one", False: "two"})
            >>> list(graph.traverse_roots())
            [0, 4, 5]
        """
        for node in self.nodes:
            if not self.parent_nodes(node):
                yield node

    def traverse_nodes_depth_first(
        self, start_node: Node
    ) -> Generator[Node, None, None]:
        """Traverses the nodes in deptch-first order that are reachable from
        the given node. In other words, the spanning graph of the given
        node is traversed in depth-first order.

        See also https://en.wikipedia.org/wiki/Depth-first_search

        Args:
            start_node: The node whose spanning graph is to be traversed
                depth first.

        Yields:
            One node of the spanning graph in depth-first order after
            another.

        Raises:
            AssertionError: If the node ``start_node`` does not exist. Note that
                the error is raised when you start to consume the returned
                generator and not when you invoke the method.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> list(graph.traverse_nodes_depth_first(0))
            [0]

            >>> list(graph.traverse_nodes_depth_first(1))
            [1]

            >>> list(graph.traverse_nodes_depth_first(2))
            [2]

            >>> list(graph.traverse_nodes_depth_first(3))
            [1, 2, 3]

            >>> list(graph.traverse_nodes_depth_first(4))
            [1, 2, 3, 4]

            >>> list(graph.traverse_nodes_depth_first(5))
            [1, 2, 5]

            >>> list(graph.traverse_nodes_depth_first(42)) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(start_node)
        nodes_to_visit = [start_node]
        reverse_depth_first_nodes = []
        while nodes_to_visit:
            current_node = nodes_to_visit.pop()
            for child_node in self.child_nodes(current_node):
                nodes_to_visit.append(child_node)
            reverse_depth_first_nodes.append(current_node)
        yielded_nodes = set()  # type: Set[Node]
        for next_node in reversed(reverse_depth_first_nodes):
            if next_node not in yielded_nodes:
                yield next_node
                yielded_nodes.add(next_node)

    def construct_spanning_graph(
        self, node: Node
    ) -> "LabeledDirectedAcyclicGraph[Node, NodeLabel, EdgeLabel]":
        """Constructs the spanning graph of the given node, that is,
        the sub-graph containing all nodes and edges reachable from the given
        node.

        Args:
            node: The node to construct the spanning graph of.

        Returns:
            The spanning graph of ``node``.

        Raises:
            AssertionError: If the node ``node`` does not exist.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> graph.add_node("zero", 0, {})
            >>> graph.add_node("one", 1, {})
            >>> graph.add_node("two", 2, {})
            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> graph.add_node("five", 5, {True: "one", False: "two"})

            >>> zero_span = graph.construct_spanning_graph(0)
            >>> (list(zero_span.nodes_with_labels), list(zero_span.edges_with_labels))
            ([(0, 'zero')], [])

            >>> one_span = graph.construct_spanning_graph(1)
            >>> (list(one_span.nodes_with_labels), list(one_span.edges_with_labels))
            ([(1, 'one')], [])

            >>> two_span = graph.construct_spanning_graph(2)
            >>> (list(two_span.nodes_with_labels), list(two_span.edges_with_labels))
            ([(2, 'two')], [])

            >>> three_span = graph.construct_spanning_graph(3)
            >>> (list(three_span.nodes_with_labels), list(three_span.edges_with_labels))
            ([(1, 'one'), (2, 'two'), (3, 'three')], [((3, 1), True), ((3, 2), False)])

            >>> four_span = graph.construct_spanning_graph(4)
            >>> (list(four_span.nodes_with_labels), list(four_span.edges_with_labels))
            ([(1, 'one'), (2, 'two'), (3, 'three'), (4, 'four')], [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False)])

            >>> five_span = graph.construct_spanning_graph(5)
            >>> (list(five_span.nodes_with_labels), list(five_span.edges_with_labels))
            ([(1, 'one'), (2, 'two'), (5, 'five')], [((5, 1), True), ((5, 2), False)])

            >>> graph.construct_spanning_graph(42) # doctest: +ELLIPSIS
            Traceback (most recent call last):
            AssertionError: The object '42' is not a node in the directed acyclic graph '<labeled_directed_acyclic_graph.LabeledDirectedAcyclicGraph object at 0x...>'
        """
        self._assert_existence_of_node(node)
        result = (
            LabeledDirectedAcyclicGraph()
        )  # type: LabeledDirectedAcyclicGraph[Node, NodeLabel, EdgeLabel]
        for current_node in self.traverse_nodes_depth_first(node):
            result.add_node(
                self.node_label(current_node),
                current_node,
                self.edge_label_to_child_node_label(current_node),
            )
        return result

    def traverse_nodes_layerwise_deepest_first(self) -> Generator[Node, None, None]:
        """Traverses the nodes starting at the deepest layer and moving up
        layer by layer, where the deepest layer consists of all leaf nodes and
        the shallowest of all root nodes.

        Note that this is not a depth-first traversal.

        Yields:
            One node after another starting at the deepest layer and moving up
            layer by layer.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            []

            >>> graph.add_node("zero", 0, {})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0]

            >>> graph.add_node("one", 1, {})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0, 1]

            >>> graph.add_node("two", 2, {})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0, 1, 2]

            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0, 1, 2, 3]

            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0, 1, 2, 3, 4]

            >>> graph.add_node("five", 5, {True: "one", False: "two"})
            >>> list(graph.traverse_nodes_layerwise_deepest_first())
            [0, 1, 2, 3, 5, 4]
        """
        visited = {node: False for node in self.nodes}
        current_layer = [node for node in self.nodes if not self.child_nodes(node)]
        while current_layer:
            next_layer = []
            for node in current_layer:
                if not visited[node]:
                    visited[node] = True
                    yield node
            # Add those unvisited parents for which all child nodes have already
            # been visited to the next layer. Unvisited parents for which this is
            # not the case will be added later as they are parents of other
            # unvisited nodes.
            for node in current_layer:
                for parent_node in self.parent_nodes(node):
                    if not visited[parent_node] and all(
                        [
                            visited[child_node]
                            for child_node in self.child_nodes(parent_node)
                        ]
                    ):
                        next_layer.append(parent_node)
            current_layer = next_layer

    def traverse_nodes_layerwise_shallowest_first(self) -> Generator[Node, None, None]:
        """Traverses the nodes starting at the shallowest layer and moving up
        layer by layer, where the shallowest layer consists of all root nodes
        and the deepest of all leaf nodes.

        Note that this is not a breadth-first traversal.

        Yields:
            One node after another starting at the shallowest layer and moving
            down layer by layer.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            []

            >>> graph.add_node("zero", 0, {})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0]

            >>> graph.add_node("one", 1, {})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0, 1]

            >>> graph.add_node("two", 2, {})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0, 1, 2]

            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0, 3, 1, 2]

            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0, 4, 3, 1, 2]

            >>> graph.add_node("five", 5, {True: "one", False: "two"})
            >>> list(graph.traverse_nodes_layerwise_shallowest_first())
            [0, 4, 5, 3, 1, 2]
        """
        visited = {node: False for node in self.nodes}
        current_layer = [node for node in self.nodes if not self.parent_nodes(node)]
        while current_layer:
            next_layer = []
            for node in current_layer:
                if not visited[node]:
                    visited[node] = True
                    yield node
            # Add those unvisited parents for which all child nodes have already
            # been visited to the next layer. Unvisited parents for which this is
            # not the case will be added later as they are parents of other
            # unvisited nodes.
            for node in current_layer:
                for child_node in self.child_nodes(node):
                    if not visited[child_node] and all(
                        [
                            visited[parent_node]
                            for parent_node in self.parent_nodes(child_node)
                        ]
                    ):
                        next_layer.append(child_node)
            current_layer = next_layer

    def traverse_forest_of_layerwise_deepest_first_spanning_graphs(
        self,
    ) -> Generator[
        Tuple[Node, "LabeledDirectedAcyclicGraph[Node, NodeLabel, EdgeLabel]"],
        None,
        None,
    ]:
        """Traverses the nodes and their spanning graphs starting at the
        deepest layer and moving up layer by layer, where the deepest layer
        consists of all leaf nodes and the shallowest of all root nodes.

        Yields:
            One pair of node and corresponding spanning graph after
            another starting at the deepest layer and moving up layer by layer.

        Examples:
            >>> graph = LabeledDirectedAcyclicGraph() # type: LabelledDirectedAcyclicGraph[int, str, bool]
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            []

            >>> graph.add_node("zero", 0, {})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], [])]

            >>> graph.add_node("one", 1, {})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], []), (1, [(1, 'one')], [])]

            >>> graph.add_node("two", 2, {})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], []), (1, [(1, 'one')], []), (2, [(2, 'two')], [])]

            >>> graph.add_node("three", 3, {True: "one", False: "two"})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], []), (1, [(1, 'one')], []), (2, [(2, 'two')], []), (3, [(1, 'one'), (2, 'two'), (3, 'three')], [((3, 1), True), ((3, 2), False)])]

            >>> graph.add_node("four", 4, {True: "one", False: "three"})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], []), (1, [(1, 'one')], []), (2, [(2, 'two')], []), (3, [(1, 'one'), (2, 'two'), (3, 'three')], [((3, 1), True), ((3, 2), False)]), (4, [(1, 'one'), (2, 'two'), (3, 'three'), (4, 'four')], [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False)])]

            >>> graph.add_node("five", 5, {True: "one", False: "two"})
            >>> [(node, list(span.nodes_with_labels), list(span.edges_with_labels)) for node, span in graph.traverse_forest_of_layerwise_deepest_first_spanning_graphs()]
            [(0, [(0, 'zero')], []), (1, [(1, 'one')], []), (2, [(2, 'two')], []), (3, [(1, 'one'), (2, 'two'), (3, 'three')], [((3, 1), True), ((3, 2), False)]), (5, [(1, 'one'), (2, 'two'), (5, 'five')], [((5, 1), True), ((5, 2), False)]), (4, [(1, 'one'), (2, 'two'), (3, 'three'), (4, 'four')], [((3, 1), True), ((3, 2), False), ((4, 1), True), ((4, 3), False)])]
        """
        for node in self.traverse_nodes_layerwise_deepest_first():
            yield (node, self.construct_spanning_graph(node))

    # private

    def _assert_existence_of_node_label(self, label: NodeLabel) -> None:
        assert label in self._label_to_node, (
            "There is no node with label '"
            + str(label)
            + "' in the graph '"
            + str(self)
            + "'"
        )

    def _assert_existence_of_node(self, node: Node) -> None:
        assert node in self._node_to_label, (
            "The object '"
            + str(node)
            + "' is not a node in the directed acyclic graph '"
            + str(self)
            + "'"
        )

    def _assert_existence_of_edge_label(self, node: Node, label: EdgeLabel) -> None:
        assert label in self._node_to_label_to_child_node[node], (
            "There is no edge from the node '"
            + str(node)
            + "' with label '"
            + str(label)
            + "' in the graph '"
            + str(self)
            + "'"
        )

    def _assert_existence_of_edge(self, edge: Tuple[Node, Node]) -> None:
        assert edge in self._edge_to_label, (
            "The object '"
            + LabeledDirectedAcyclicGraph.format_edge(edge)
            + "' is not an edge in the directed acyclic graph '"
            + str(self)
            + "'"
        )
