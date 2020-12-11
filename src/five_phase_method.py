"""Implementation of the five-phase method based on Radiance for the
computation of annual illuminance and images."""

# pylint: disable=invalid-name

from labeled_directed_acyclic_graph import LabeledDirectedAcyclicGraph
from data_flow_graph_executor import DataFlowGraphExecutor, Transformation
from data_flow_graph import (
    DataFlowLabel,
    Input,
    DataFormat,
    Output,
    TransmissionKind,
    ComputationLabel,
    SaveToFile,
    File,
    DataFlowGraph,
    Series,
    ShellCommand,
    Computation,
)
from pathlib import Path


# Ulimit needs to be set to a high enough value so that multiple files can be
# opened simulataneously.
# For the present simulation, maximum number of open files will be 5165.
# This command is applicable to Unix-like systems only. On Windows the number of
# simultaneous files is limited to 512.
# ulimit -n 9999
def compute(
    *,
    annual_illuminance_output_path: Path,
    bsdf_xml_path: Path,
    glazing_aperture_geometry_material_name: str,
    images_output_path: Path,
    materialless_glazing_aperture_geometry_path: Path,
    outside_scene_path: Path,
    points_path: Path,
    room_material_path: Path,
    room_scene_path: Path,
    tensor_tree_path: Path,
    typical_meteorological_year_weather_path: Path,
    view_path: Path
) -> None:
    """Computes the annual illuminance and the images with the five-phase
    method for the given inputs."""
    data_flow_graph = FivePhaseMethod.graph(
        annual_illuminance_output_path=annual_illuminance_output_path,
        bsdf_xml_path=bsdf_xml_path,
        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
        images_output_path=images_output_path,
        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
        outside_scene_path=outside_scene_path,
        points_path=points_path,
        room_material_path=room_material_path,
        room_scene_path=room_scene_path,
        tensor_tree_path=tensor_tree_path,
        typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
        view_path=view_path,
    )
    DataFlowGraphExecutor(
        data_flow_graph, transformation=Transformation.SPLITS
    ).execute()


def compute_annual_illuminance(
    *,
    annual_illuminance_output_path: Path,
    bsdf_xml_path: Path,
    glazing_aperture_geometry_material_name: str,
    materialless_glazing_aperture_geometry_path: Path,
    outside_scene_path: Path,
    points_path: Path,
    room_material_path: Path,
    room_scene_path: Path,
    tensor_tree_path: Path,
    typical_meteorological_year_weather_path: Path
) -> None:
    """Computes the annual illuminance with the five-phase method for the given
    inputs."""
    data_flow_graph = FivePhaseMethod.save_annual_illuminance_graph(
        annual_illuminance_output_path=annual_illuminance_output_path,
        bsdf_xml_path=bsdf_xml_path,
        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
        outside_scene_path=outside_scene_path,
        points_path=points_path,
        room_material_path=room_material_path,
        room_scene_path=room_scene_path,
        tensor_tree_path=tensor_tree_path,
        typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
    )
    DataFlowGraphExecutor(
        data_flow_graph, transformation=Transformation.DIAMONDS
    ).execute()


def compute_images(
    *,
    bsdf_xml_path: Path,
    glazing_aperture_geometry_material_name: str,
    images_output_path: Path,
    materialless_glazing_aperture_geometry_path: Path,
    outside_scene_path: Path,
    room_material_path: Path,
    room_scene_path: Path,
    tensor_tree_path: Path,
    typical_meteorological_year_weather_path: Path,
    view_path: Path
) -> None:
    """Computes the images with the five-phase method for the given inputs."""
    data_flow_graph = FivePhaseMethod.save_images_graph(
        bsdf_xml_path=bsdf_xml_path,
        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
        images_output_path=images_output_path,
        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
        outside_scene_path=outside_scene_path,
        room_material_path=room_material_path,
        room_scene_path=room_scene_path,
        tensor_tree_path=tensor_tree_path,
        typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
        view_path=view_path,
    )
    DataFlowGraphExecutor(
        data_flow_graph, transformation=Transformation.DIAMONDS
    ).execute()


def _frozen_octree_from_material_and_scene() -> Computation:
    return ShellCommand(
        formatted_command="oconv -f {material} {scene}",
        label_to_input={
            DataFlowLabel("material"): Input(data_format=DataFormat.MAT),
            DataFlowLabel("scene"): Input(data_format=DataFormat.RAD),
        },
        output=Output(data_format=DataFormat.OCT),
    )


def _frozen_octree_from_material_and_scene_and_geometry() -> Computation:
    return ShellCommand(
        formatted_command="oconv -f {material} {scene} {geometry}",
        label_to_input={
            DataFlowLabel("material"): Input(data_format=DataFormat.MAT),
            DataFlowLabel("scene"): Input(data_format=DataFormat.RAD),
            DataFlowLabel("geometry"): Input(data_format=DataFormat.RAD),
        },
        output=Output(data_format=DataFormat.OCT),
    )


class Root:
    # 5.2 Annual sky-matrix using TMY weather data

    ##Create sky-vectors
    ##Annual sky-matrix
    # epw2wea assets/USA_NY_New.York-Central.Park.725033_TMY3m.epw assets/NYC.wea
    # aka `assets/NYC.wea`
    @staticmethod
    def typical_meteorological_year_weather_graph(
        *, typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("typical_meteorological_year_weather",)),
            File(typical_meteorological_year_weather_path),
            {},
        )

    ###Create an annual sky matrix with 145 patches.
    # gendaymtx -m 1 assets/NYC.wea > skyVectors/NYC.smx
    # aka `skyVectors/NYC.smx`
    # Something similar is generated by `genMtx#sky`.
    @staticmethod
    def annual_sky_matrix_graph(
        *, typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("annual_sky_matrix",)),
            ShellCommand(
                formatted_command="gendaymtx -m 1 {weather}",
                label_to_input={
                    DataFlowLabel("weather"): Input(data_format=DataFormat.WEA)
                },
                output=Output(data_format=DataFormat.SMX),
            ),
            {
                DataFlowLabel(
                    "weather"
                ): Root.typical_meteorological_year_weather_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                )
            },
        )

    # 7.1 Flux Transfer Matrices
    # 7.1.1 View Matrix

    # aka `materials.rad`
    @staticmethod
    def room_material_graph(*, room_material_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("room_material",)), File(room_material_path), {}
        )

    # aka `room.rad`
    @staticmethod
    def room_scene_graph(*, room_scene_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("room_scene",)), File(room_scene_path), {}
        )

    # aka `points.txt`
    @staticmethod
    def points_graph(*, points_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("points",)), File(points_path), {}
        )

    # glazing aperture geometry
    # aka `objects/GlazingVmtx.rad`
    @staticmethod
    def materialless_glazing_aperture_geometry_graph(
        *, materialless_glazing_aperture_geometry_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("materialless_glazing_aperture_geometry",)),
            File(materialless_glazing_aperture_geometry_path),
            {},
        )

    @staticmethod
    def glazing_aperture_geometry_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("glazing_aperture_geometry",)),
            ShellCommand(
                formatted_command=(
                    "echo -e 'void glow {material_name} 0 0 4 1 1 1 0\n'".format(
                        material_name=glazing_aperture_geometry_material_name
                    )
                    + " && cat {materialless_geometry}"
                ),
                label_to_input={
                    DataFlowLabel("materialless_geometry"): Input(
                        data_format=DataFormat.RAD
                    )
                },
                output=Output(data_format=DataFormat.RAD),
            ),
            {
                DataFlowLabel(
                    "materialless_geometry"
                ): Root.materialless_glazing_aperture_geometry_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path
                )
            },
        )

    # aka `views/south.vf`
    @staticmethod
    def view_graph(*, view_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(ComputationLabel(("view",)), File(view_path), {})

    # aka `matrices/tmtx/blinds.xml` or `blinds/blinds.xml`
    @staticmethod
    def bsdf_xml_graph(*, bsdf_xml_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("bsdf_xml",)), File(bsdf_xml_path), {}
        )

    @staticmethod
    def tensor_tree_graph(*, tensor_tree_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("tensor_tree",)), File(tensor_tree_path), {}
        )

    @staticmethod
    def func_cal_graph() -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("func_cal",)),
            ShellCommand(
                formatted_command=(
                    "echo -e '" " ux = 0;\n" " uy = 0;\n" " uz = 1;\n" " thick = 0;" "'"
                ),
                label_to_input={},
                output=Output(data_format=DataFormat.CAL),
            ),
            {},
        )

    # aka `blinds/blindsWithProxy.rad`
    @staticmethod
    def bsdf_with_proxy_rad_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        tensor_tree_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("bsdf_with_proxy_rad",)),
            ShellCommand(
                formatted_command=(
                    "echo -e '"
                    + " void BSDF {material_name}\n".format(
                        material_name=glazing_aperture_geometry_material_name
                    )
                    + " 6 thick {tensor_tree} ux uy uz {func_cal}\n"
                    " 0\n"
                    " 0\n"
                    "'"
                    " && cat {materialless_geometry}"
                ),
                label_to_input={
                    DataFlowLabel("tensor_tree"): Input(data_format=DataFormat.XML),
                    DataFlowLabel("func_cal"): Input(
                        data_format=DataFormat.CAL,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("materialless_geometry"): Input(
                        data_format=DataFormat.RAD
                    ),
                },
                output=Output(data_format=DataFormat.RAD),
            ),
            {
                DataFlowLabel("tensor_tree"): Root.tensor_tree_graph(
                    tensor_tree_path=tensor_tree_path
                ),
                DataFlowLabel("func_cal"): Root.func_cal_graph(),
                DataFlowLabel(
                    "materialless_geometry"
                ): Root.materialless_glazing_aperture_geometry_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path
                ),
            },
        )

    # 7.1.3 Daylight Matrix

    # sky and ground definition
    # aka `skyDomes/skyglow.rad`
    @staticmethod
    def outside_scene_graph(*, outside_scene_path: Path) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("outside_scene",)), File(outside_scene_path), {}
        )

    # https://github.com/NREL/Radiance/blob/master/src/cal/cal/reinsrc.cal
    # aka `reinsrc.cal`
    @staticmethod
    def reinhart_sky_directions_calculation_graph() -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("reinhart_sky_directions_calculation",)),
            File(Path("calculations/reinsrc.cal")),
            {},
        )

    # https://github.com/NREL/Radiance/blob/master/src/cal/cal/reinhart.cal
    # aka `reinhart.cal`
    @staticmethod
    def reinhart_high_density_sky_patches_calculation_graph() -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("reinhart_high_density_sky_patches_calculation",)),
            File(Path("calculations/reinhart.cal")),
            {},
        )

    # 5.1.3 Creating the sky-vectors

    # gendaylit 3 20 10:30EDT -m 75 -o 73.96 -a 40.78 -W 706 162 > skies/NYC_Per_DH.sky
    # aka `skies/NYC_Per_DH.sky`
    @staticmethod
    def perez_sky_graph(
        *,
        month: int,
        day: int,
        hour: float,
        time_zone: float,
        latitude: float,
        longitude: float,
        direct_normal_irradiance: float,
        diffuse_horizontal_irradiance: float
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("perez_sky",)),
            ShellCommand(
                formatted_command=(
                    "gendaylit {month:d} {day:d} {hour:f}"
                    " -a {latitude:1.2f} -o {longitude:1.2f}"
                    " -m {time_zone:1.2f}"
                    " -W {direct_normal_irradiance:1.2f} {diffuse_horizontal_irradiance:1.2f}"
                ).format(
                    month=month,
                    day=day,
                    hour=hour,
                    latitude=latitude,
                    longitude=longitude,
                    time_zone=time_zone,
                    direct_normal_irradiance=direct_normal_irradiance,
                    diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                ),
                label_to_input={},
                output=Output(data_format=DataFormat.SKY),
            ),
            {},
        )

    # genskyvec -m 1 < skies/NYC_Per_DH.sky > skyVectors/NYC_Per.vec
    # aka `skyVectors/NYC_Per.vec`
    @staticmethod
    def perez_sky_vector_graph(
        *,
        month: int,
        day: int,
        hour: float,
        time_zone: float,
        latitude: float,
        longitude: float,
        direct_normal_irradiance: float,
        diffuse_horizontal_irradiance: float
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(("perez_sky_vector",)),
            ShellCommand(
                formatted_command="genskyvec -m 1 < {sky}",
                label_to_input={
                    DataFlowLabel("sky"): Input(data_format=DataFormat.SKY)
                },
                output=Output(data_format=DataFormat.VEC),
            ),
            {
                DataFlowLabel("sky"): Root.perez_sky_graph(
                    month=month,
                    day=day,
                    hour=hour,
                    time_zone=time_zone,
                    latitude=latitude,
                    longitude=longitude,
                    direct_normal_irradiance=direct_normal_irradiance,
                    diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                )
            },
        )

    class Black:

        # aka `materialBlack.rad`
        @staticmethod
        def black_material_graph() -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(("black_material",)),
                ShellCommand(
                    formatted_command="echo -e 'void plastic black 0 0 5 0 0 0 0 0'",
                    label_to_input={},
                    output=Output(data_format=DataFormat.MAT),
                ),
                {},
            )

        # aka `roomBlack.rad`
        @staticmethod
        def black_room_scene_graph(*, room_scene_path: Path) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(("black_room_scene",)),
                ShellCommand(
                    formatted_command="sed -e 's/!xform/!xform -m black/g' {scene}",
                    label_to_input={
                        DataFlowLabel("scene"): Input(data_format=DataFormat.RAD)
                    },
                    output=Output(data_format=DataFormat.RAD),
                ),
                {
                    DataFlowLabel("scene"): Root.room_scene_graph(
                        room_scene_path=room_scene_path
                    )
                },
            )

        # Create a black octree for direct calculation.
        # oconv -f materialBlack.rad roomBlack.rad > octrees/room3phDirect.oct
        # aka `octrees/room3phDirect.oct`
        @staticmethod
        def black_room_octree_graph(*, room_scene_path: Path) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(("black_room_octree",)),
                _frozen_octree_from_material_and_scene(),
                {
                    DataFlowLabel("material"): Root.Black.black_material_graph(),
                    DataFlowLabel("scene"): Root.Black.black_room_scene_graph(
                        room_scene_path=room_scene_path
                    ),
                },
            )

        # Create an octree with opaque glazing for material map.
        # oconv -f room.mat room.rad objects/GlazingVmtx.rad > octrees/materialMap3PhaseDirect.oct
        # aka `octrees/materialMap3PhaseDirect.oct`
        @staticmethod
        def opaque_glazing_octree_graph(
            *,
            glazing_aperture_geometry_material_name: str,
            materialless_glazing_aperture_geometry_path: Path,
            room_scene_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(("opaque_glazing_octree",)),
                _frozen_octree_from_material_and_scene_and_geometry(),
                {
                    DataFlowLabel("material"): Root.Black.black_material_graph(),
                    DataFlowLabel("scene"): Root.Black.black_room_scene_graph(
                        room_scene_path=room_scene_path
                    ),
                    DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    ),
                },
            )


class ThreePhaseMethod:
    GROUP_PATH = ("three_phase_method",)

    class IlluminanceAndImage:

        # 7.2 Generating Results

        # dctimestep matrices/vmtx/v.mtx matrices/tmtx/blinds.xml matrices/dmtx/daylight.dmx skyVectors/NYC_Per.vec
        #  | rmtxop -fa -c 47.4 119.9 11.6 -
        #  > results/3ph/3ph.ill
        # aka `results/3ph/3ph.ill`
        @staticmethod
        def illuminance_graph(
            *,
            bsdf_xml_path: Path,
            glazing_aperture_geometry_material_name: str,
            outside_scene_path: Path,
            materialless_glazing_aperture_geometry_path: Path,
            points_path: Path,
            room_material_path: Path,
            room_scene_path: Path,
            month: int,
            day: int,
            hour: float,
            time_zone: float,
            latitude: float,
            longitude: float,
            direct_normal_irradiance: float,
            diffuse_horizontal_irradiance: float
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("illuminance",)),
                ShellCommand(
                    formatted_command=(
                        "dctimestep {view_matrix} {bsdf} {transmission_matrix} {sky_vector}"
                        " | rmtxop -fa -c 47.4 119.9 11.6 -"
                    ),
                    label_to_input={
                        DataFlowLabel("view_matrix"): Input(data_format=DataFormat.MTX),
                        DataFlowLabel("bsdf"): Input(
                            transmission_kind=TransmissionKind.FILE,
                            data_format=DataFormat.XML,
                        ),  # `bsdf_xml` cannot be piped for some reason. If we try do so, `dctimestep` produces the error `fatal - unexpected EOF in header`.
                        DataFlowLabel("transmission_matrix"): Input(
                            data_format=DataFormat.DMX
                        ),
                        DataFlowLabel("sky_vector"): Input(data_format=DataFormat.VEC),
                    },
                    output=Output(data_format=DataFormat.ILL),
                ),
                {
                    DataFlowLabel(
                        "view_matrix"
                    ): ThreePhaseMethod.Simulation.illuminance_view_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        points_path=points_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                    DataFlowLabel("bsdf"): Root.bsdf_xml_graph(
                        bsdf_xml_path=bsdf_xml_path
                    ),
                    DataFlowLabel(
                        "transmission_matrix"
                    ): ThreePhaseMethod.Simulation.daylight_transmission_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                    DataFlowLabel("sky_vector"): Root.perez_sky_vector_graph(
                        month=month,
                        day=day,
                        hour=hour,
                        time_zone=time_zone,
                        latitude=latitude,
                        longitude=longitude,
                        direct_normal_irradiance=direct_normal_irradiance,
                        diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                    ),
                },
            )

        @staticmethod
        def save_illuminance_graph(
            *,
            bsdf_xml_path: Path,
            glazing_aperture_geometry_material_name: str,
            materialless_glazing_aperture_geometry_path: Path,
            illuminance_output_path: Path,
            outside_scene_path: Path,
            points_path: Path,
            room_material_path: Path,
            room_scene_path: Path,
            month: int,
            day: int,
            hour: float,
            time_zone: float,
            latitude: float,
            longitude: float,
            direct_normal_irradiance: float,
            diffuse_horizontal_irradiance: float
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("save_illuminance",)),
                SaveToFile(
                    label_to_input_with_path={
                        DataFlowLabel("data"): (
                            Input(data_format=DataFormat.ILL),
                            illuminance_output_path,
                        )
                    }
                ),
                {
                    DataFlowLabel(
                        "data"
                    ): ThreePhaseMethod.IlluminanceAndImage.illuminance_graph(
                        bsdf_xml_path=bsdf_xml_path,
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        points_path=points_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                        month=month,
                        day=day,
                        hour=hour,
                        time_zone=time_zone,
                        latitude=latitude,
                        longitude=longitude,
                        direct_normal_irradiance=direct_normal_irradiance,
                        diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                    )
                },
            )

        # dctimestep -h matrices/vmtx/hdr/south%03d.hdr matrices/tmtx/blinds.xml
        # matrices/dmtx/daylight.dmx skyVectors/NYC_Per.vec > results/3ph/3ph.hdr
        # aka `results/3ph/3ph.hdr`
        @staticmethod
        def image_graph(
            *,
            bsdf_xml_path: Path,
            glazing_aperture_geometry_material_name: str,
            materialless_glazing_aperture_geometry_path: Path,
            outside_scene_path: Path,
            room_material_path: Path,
            view_path: Path,
            room_scene_path: Path,
            month: int,
            day: int,
            hour: float,
            time_zone: float,
            latitude: float,
            longitude: float,
            direct_normal_irradiance: float,
            diffuse_horizontal_irradiance: float
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("image",)),
                ShellCommand(
                    formatted_command=(
                        "dctimestep -h {view_matrix} {bsdf} {transmission_matrix} {sky_vector}"
                    ),
                    label_to_input={
                        DataFlowLabel("view_matrix"): Input(
                            data_format=DataFormat.HDR,
                            transmission_kind=TransmissionKind.FILE,
                            series=Series(format="%03d"),
                        ),
                        DataFlowLabel("bsdf"): Input(data_format=DataFormat.XML),
                        DataFlowLabel("transmission_matrix"): Input(
                            data_format=DataFormat.DMX
                        ),
                        DataFlowLabel("sky_vector"): Input(data_format=DataFormat.VEC),
                    },
                    output=Output(data_format=DataFormat.HDR),
                ),
                {
                    DataFlowLabel(
                        "view_matrix"
                    ): ThreePhaseMethod.Simulation.images_view_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                        view_path=view_path,
                    ),
                    DataFlowLabel("bsdf"): Root.bsdf_xml_graph(
                        bsdf_xml_path=bsdf_xml_path
                    ),
                    DataFlowLabel(
                        "transmission_matrix"
                    ): ThreePhaseMethod.Simulation.daylight_transmission_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        room_scene_path=room_scene_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_material_path=room_material_path,
                    ),
                    DataFlowLabel("sky_vector"): Root.perez_sky_vector_graph(
                        month=month,
                        day=day,
                        hour=hour,
                        time_zone=time_zone,
                        latitude=latitude,
                        longitude=longitude,
                        direct_normal_irradiance=direct_normal_irradiance,
                        diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                    ),
                },
            )

        @staticmethod
        def save_image_graph(
            *,
            bsdf_xml_path: Path,
            glazing_aperture_geometry_material_name: str,
            image_output_path: Path,
            outside_scene_path: Path,
            room_material_path: Path,
            view_path: Path,
            materialless_glazing_aperture_geometry_path: Path,
            room_scene_path: Path,
            month: int,
            day: int,
            hour: float,
            time_zone: float,
            latitude: float,
            longitude: float,
            direct_normal_irradiance: float,
            diffuse_horizontal_irradiance: float
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("save_image",)),
                SaveToFile(
                    label_to_input_with_path={
                        DataFlowLabel("data"): (
                            Input(data_format=DataFormat.HDR),
                            image_output_path,
                        )
                    }
                ),
                {
                    DataFlowLabel(
                        "data"
                    ): ThreePhaseMethod.IlluminanceAndImage.image_graph(
                        bsdf_xml_path=bsdf_xml_path,
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        room_material_path=room_material_path,
                        view_path=view_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_scene_path=room_scene_path,
                        month=month,
                        day=day,
                        hour=hour,
                        time_zone=time_zone,
                        latitude=latitude,
                        longitude=longitude,
                        direct_normal_irradiance=direct_normal_irradiance,
                        diffuse_horizontal_irradiance=diffuse_horizontal_irradiance,
                    )
                },
            )

    class Simulation:
        # PART1: THREE-PHASE METHOD SIMULATION.
        ##Create octree
        # oconv -f room.mat room.rad > octrees/room3ph.oct
        # aka `octrees/room3ph.oct`
        @staticmethod
        def room_octree_graph(
            *, room_material_path: Path, room_scene_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("room_octree",)),
                _frozen_octree_from_material_and_scene(),
                {
                    DataFlowLabel("material"): Root.room_material_graph(
                        room_material_path=room_material_path
                    ),
                    DataFlowLabel("scene"): Root.room_scene_graph(
                        room_scene_path=room_scene_path
                    ),
                },
            )

        # 7.1.1.1 View Matrix for Illuminance

        ##View matrix for Illuminance (The value for -y 100 is derived from the 100 grid points inside the file points.txt).
        # rfluxmtx -v -I+ -ab 4 -ad 5000 -lw 0.0002 -n 16 -y 100 - objects/GlazingVmtx.rad -i octrees/room3ph.oct < points.txt > matrices/vmtx/v.mtx
        # aka `matrices/vmtx/v.mtx`
        @staticmethod
        def illuminance_view_matrix_graph(
            *,
            glazing_aperture_geometry_material_name: str,
            materialless_glazing_aperture_geometry_path: Path,
            points_path: Path,
            room_material_path: Path,
            room_scene_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(
                    ThreePhaseMethod.GROUP_PATH + ("illuminance_view_matrix",)
                ),
                ShellCommand(
                    formatted_command="rfluxmtx -v -I+ -ab 4 -ad 5000 -lw 0.0002 -n 16 -y 100 - {geometry} -i {octree} < {points}",
                    label_to_input={
                        DataFlowLabel("geometry"): Input(
                            data_format=DataFormat.RAD, repetition_count=2
                        ),
                        DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                        DataFlowLabel("points"): Input(data_format=DataFormat.PTS),
                    },
                    output=Output(data_format=DataFormat.MTX),
                ),
                {
                    DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    ),
                    DataFlowLabel(
                        "octree"
                    ): ThreePhaseMethod.Simulation.room_octree_graph(
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                    DataFlowLabel("points"): Root.points_graph(points_path=points_path),
                },
            )

        # 7.1.1.2 View Matrix for Image-based simulations

        ##View matrix for Images.
        # vwrays -vf views/south.vf -x 400 -y 400 -pj 0.7 -c 9 -ff
        #  | rfluxmtx -v -ffc `vwrays -vf views/south.vf -x 400 -y 400 -d`
        #  -o matrices/vmtx/hdr/south%03d.hdr -ab 4 -ad 1000 -lw 1e-4 -c 9 -n 16
        #  - objects/GlazingVmtx.rad -i octrees/room3ph.oct
        # aka `matrices/vmtx/hdr/south%03d.hdr`
        @staticmethod
        def images_view_matrix_graph(
            *,
            glazing_aperture_geometry_material_name: str,
            materialless_glazing_aperture_geometry_path: Path,
            room_material_path: Path,
            room_scene_path: Path,
            view_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("images_view_matrix",)),
                ShellCommand(
                    formatted_command=(
                        "vwrays -vf {view} -x 400 -y 400 -pj 0.7 -c 9 -ff"
                        " | rfluxmtx -v -ffc `vwrays -vf {view} -x 400 -y 400 -d`"
                        " -o {output} -ab 4 -ad 1000 -lw 1e-4 -c 9 -n 16"
                        " - {geometry} -i {octree}"
                    ),
                    label_to_input={
                        # TODO Check whether the inputs may be streamed. Note that `view` is used twice in the command.
                        DataFlowLabel("view"): Input(
                            data_format=DataFormat.VF,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                        DataFlowLabel("geometry"): Input(
                            data_format=DataFormat.RAD,
                            # repetition_count=2,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                        DataFlowLabel("octree"): Input(
                            data_format=DataFormat.OCT,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                    },
                    output=Output(
                        data_format=DataFormat.HDR, series=Series(format="%03d")
                    ),
                ),
                {
                    DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                    DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    ),
                    DataFlowLabel(
                        "octree"
                    ): ThreePhaseMethod.Simulation.room_octree_graph(
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                },
            )

        # 7.1.2 Transmission Matrix

        ##D matrix
        # rfluxmtx -v -ff -ab 4 -ad 1000 -lw 0.001 -c 1000 -n 16 objects/GlazingVmtx.rad
        # skyDomes/skyglow.rad -i octrees/room3ph.oct > matrices/dmtx/daylight.dmx
        # aka `matrices/dmtx/daylight.dmx`
        @staticmethod
        def daylight_transmission_matrix_graph(
            *,
            glazing_aperture_geometry_material_name: str,
            outside_scene_path: Path,
            materialless_glazing_aperture_geometry_path: Path,
            room_material_path: Path,
            room_scene_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(
                    ThreePhaseMethod.GROUP_PATH + ("daylight_transmission_matrix",)
                ),
                ShellCommand(
                    formatted_command="rfluxmtx -v -ff -ab 4 -ad 1000 -lw 0.001 -c 1000 -n 16 {geometry} {scene} -i {octree}",
                    label_to_input={
                        DataFlowLabel("geometry"): Input(data_format=DataFormat.RAD),
                        DataFlowLabel("scene"): Input(
                            data_format=DataFormat.RAD, repetition_count=2
                        ),
                        DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                    },
                    output=Output(data_format=DataFormat.DMX),
                ),
                {
                    DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    ),
                    DataFlowLabel("scene"): Root.outside_scene_graph(
                        outside_scene_path=outside_scene_path
                    ),
                    DataFlowLabel(
                        "octree"
                    ): ThreePhaseMethod.Simulation.room_octree_graph(
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                },
            )

        # 7.2 Generating Results

        ###Illuminance
        # dctimestep matrices/vmtx/v.mtx matrices/tmtx/blinds.xml
        # matrices/dmtx/daylight.dmx skyVectors/NYC.smx
        #  | rmtxop -fa -t -c 47.4 119.9 11.6 -
        #  > results/5ph/3ph/3phAnnual.ill
        # aka `results/3ph/3phAnnual.ill`
        # annual_illuminance = File(
        #     ComputationLabel(GROUP_PATH + ("annual_illuminance",)), path="backup/annual_illuminance.ill"
        # )
        # data_flow_graph.add_node(annual_illuminance)
        @staticmethod
        def annual_illuminance_graph(
            *,
            glazing_aperture_geometry_material_name: str,
            bsdf_xml_path: Path,
            outside_scene_path: Path,
            materialless_glazing_aperture_geometry_path: Path,
            points_path: Path,
            room_material_path: Path,
            room_scene_path: Path,
            typical_meteorological_year_weather_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("annual_illuminance",)),
                ShellCommand(
                    formatted_command=(
                        "dctimestep {view_matrix} {bsdf} {transmission_matrix} {sky_matrix}"
                        " | rmtxop -fa -t -c 47.4 119.9 11.6 -"
                    ),
                    label_to_input={
                        DataFlowLabel("view_matrix"): Input(data_format=DataFormat.MTX),
                        DataFlowLabel("bsdf"): Input(
                            transmission_kind=TransmissionKind.FILE,
                            data_format=DataFormat.XML,
                        ),  # TODO `bsdf_xml` cannot be piped for some reason. If we try do so, `dctimestep` produces the error `fatal - unexpected EOF in header`. Is this also the case for other usages of `dctimestep`?
                        DataFlowLabel("transmission_matrix"): Input(
                            data_format=DataFormat.DMX
                        ),
                        DataFlowLabel("sky_matrix"): Input(data_format=DataFormat.SMX),
                    },
                    output=Output(data_format=DataFormat.ILL),
                ),
                {
                    DataFlowLabel(
                        "view_matrix"
                    ): ThreePhaseMethod.Simulation.illuminance_view_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        points_path=points_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                    DataFlowLabel("bsdf"): Root.bsdf_xml_graph(
                        bsdf_xml_path=bsdf_xml_path
                    ),
                    DataFlowLabel(
                        "transmission_matrix"
                    ): ThreePhaseMethod.Simulation.daylight_transmission_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                    ),
                    DataFlowLabel("sky_matrix"): Root.annual_sky_matrix_graph(
                        typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                    ),
                },
            )

        ###Results for views/south.vf
        ## Command from appendix E
        # dctimestep -o results/5ph/3ph/hdr/south%04d.hdr matrices/vmtx/hdr/south%03d.hdr
        # matrices/tmtx/blinds.xml matrices/dmtx/daylight.dmx skyVectors/NYC.smx
        # aka `results/3ph/hdr/south%04d.hdr`
        # images = File(
        #     ComputationLabel(GROUP_PATH + ("images",)),
        #     path="backup/tmp/tmpvwnq6e92/file0_of_south_3ph_hdr_%04d.hdr",
        # )
        # data_flow_graph.add_node(images)
        @staticmethod
        def images_graph(
            *,
            bsdf_xml_path: Path,
            glazing_aperture_geometry_material_name: str,
            outside_scene_path: Path,
            materialless_glazing_aperture_geometry_path: Path,
            room_material_path: Path,
            room_scene_path: Path,
            typical_meteorological_year_weather_path: Path,
            view_path: Path
        ) -> DataFlowGraph:
            return DataFlowGraph.construct(
                ComputationLabel(ThreePhaseMethod.GROUP_PATH + ("images",)),
                ShellCommand(
                    formatted_command=(
                        "dctimestep -o {output} {view_matrix} {bsdf} {transmission_matrix} {sky_matrix}"
                    ),
                    label_to_input={
                        DataFlowLabel("view_matrix"): Input(
                            data_format=DataFormat.HDR,
                            transmission_kind=TransmissionKind.FILE,
                            series=Series(format="%03d"),
                        ),
                        # TODO Check whether the following inputs may be streamed
                        DataFlowLabel("bsdf"): Input(
                            data_format=DataFormat.XML,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                        DataFlowLabel("transmission_matrix"): Input(
                            data_format=DataFormat.DMX,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                        DataFlowLabel("sky_matrix"): Input(
                            data_format=DataFormat.SMX,
                            transmission_kind=TransmissionKind.FILE,
                        ),
                    },
                    output=Output(
                        data_format=DataFormat.HDR, series=Series(format="%04d")
                    ),
                ),
                {
                    DataFlowLabel(
                        "view_matrix"
                    ): ThreePhaseMethod.Simulation.images_view_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_material_path=room_material_path,
                        room_scene_path=room_scene_path,
                        view_path=view_path,
                    ),
                    DataFlowLabel("bsdf"): Root.bsdf_xml_graph(
                        bsdf_xml_path=bsdf_xml_path
                    ),
                    DataFlowLabel(
                        "transmission_matrix"
                    ): ThreePhaseMethod.Simulation.daylight_transmission_matrix_graph(
                        glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                        outside_scene_path=outside_scene_path,
                        materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                        room_scene_path=room_scene_path,
                        room_material_path=room_material_path,
                    ),
                    DataFlowLabel("sky_matrix"): Root.annual_sky_matrix_graph(
                        typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                    ),
                },
            )

        # 7.3 Parametric simulations by reusing phases

        # TODO Are the following calculations needed somewhere?
        # dctimestep matrices/vmtx/v.mtx matrices/tmtx/ven0.xml matrices/dmtx/daylight.dmx
        # skyVectors/NYC_Per.vec | rmtxop -fa -c 47.4 119.9 11.6 - >
        # results/3ph/3phVen0.ill

        # dctimestep matrices/vmtx/hdr/south%03d.hdr matrices/tmtx/ven0.xml
        # matrices/dmtx/daylight.dmx skyVectors/NYC_Per.vec > results/3ph/3phVen0.hdr

        # 7.4 Simulating a vertically adjustable shading system

        # 7.5 Simulating non-coplanar shading systems

        # Appendix E. Commands for the Five-Phase and Six-Phase Methods
        # E.1 Commands for the Five-Phase Method


class DirectSolarPartOfThreePhaseMethodSimulation:
    GROUP_PATH = ("direct_solar_part",)

    # PART2: DIRECT SOLAR PART OF THE THREE-PHASE SIMULATION.

    # Direct View matrix for Illuminance
    # rfluxmtx -v -I+ -ab 1 -ad 5000 -lw 0.0002 -n 16 -y 100 - objects/GlazingVmtx.rad -i octrees/room3phDirect.oct < points.txt > matrices/vmtxd/v.mtx
    # aka `matrices/vmtxd/v.mtx`
    @staticmethod
    def illuminance_view_matrix_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        points_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("illuminance_view_matrix",)
            ),
            ShellCommand(
                formatted_command="rfluxmtx -v -I+ -ab 1 -ad 5000 -lw 0.0002 -n 16 -y 100 - {geometry} -i {octree} < {points}",
                label_to_input={
                    DataFlowLabel("geometry"): Input(
                        data_format=DataFormat.RAD, repetition_count=2
                    ),
                    DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                    DataFlowLabel("points"): Input(data_format=DataFormat.PTS),
                },
                output=Output(data_format=DataFormat.MTX),
            ),
            {
                DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                ),
                DataFlowLabel("octree"): Root.Black.black_room_octree_graph(
                    room_scene_path=room_scene_path
                ),
                DataFlowLabel("points"): Root.points_graph(points_path=points_path),
            },
        )

    ##Direct View matrix for Images.
    # vwrays -vf views/south.vf -x 400 -y 400 -pj 0.7 -c 9 -ff
    #  | rfluxmtx -v -ffc -i `vwrays -vf views/south.vf -x 400 -y 400 -d`
    #  -o matrices/vmtxd/hdrIllum/south%03d.hdr -ab 1 -ad 1000 -lw 1e-4 -c 9 -n 16
    #  - objects/GlazingVmtx.rad -i octrees/room3phDirect.oct
    # aka `matrices/vmtxd/hdrIllum/south%03d.hdr`
    @staticmethod
    def images_view_matrix_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        room_scene_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("images_view_matrix",)
            ),
            ShellCommand(
                formatted_command=(
                    "vwrays -vf {view} -x 400 -y 400 -pj 0.7 -c 9 -ff"
                    " | rfluxmtx -v -ffc -i"
                    " `vwrays -vf {view} -x 400 -y 400 -d`"
                    " -o {output} -ab 1 -ad 1000 -lw 1e-4 -c 9 -n 16"
                    " - {geometry} -i {octree}"
                ),
                label_to_input={
                    # TODO Check whether the inputs `view` and `octree` may be streamed. Note that `view` occurs twice in the command!
                    DataFlowLabel("view"): Input(
                        data_format=DataFormat.VF,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("geometry"): Input(
                        data_format=DataFormat.RAD,
                        # repetition_count=2,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("octree"): Input(
                        data_format=DataFormat.OCT,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                },
                output=Output(data_format=DataFormat.HDR, series=Series(format="%03d")),
            ),
            {
                DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                ),
                DataFlowLabel("octree"): Root.Black.black_room_octree_graph(
                    room_scene_path=room_scene_path
                ),
            },
        )

    # rpict -x 400 -y 400 -ps 1 -av 0.31831 0.31831 0.31831 -ab 0 -vf views/south.vf octrees/materialMap3PhaseDirect.oct > matrices/vmtxd/materialMapSouth.hdr
    # aka `matrices/vmtxd/materialMapSouth.hdr`
    @staticmethod
    def material_map_image_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("material_map_image",)
            ),
            ShellCommand(
                formatted_command="rpict -x 400 -y 400 -ps 1 -av 0.31831 0.31831 0.31831 -ab 0 -vf {view} {octree}",
                label_to_input={
                    DataFlowLabel("view"): Input(data_format=DataFormat.VF),
                    DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                },
                output=Output(data_format=DataFormat.HDR),
            ),
            {
                DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                DataFlowLabel("octree"): Root.Black.opaque_glazing_octree_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                ),
            },
        )

    # Calculate luminance-based images from illuminance based images using material
    # maps.
    # for idx in {00..145}
    # do
    #     pcomb -h -e 'ro=ri(1)*ri(2);go=gi(1)*gi(2);bo=bi(1)*bi(2)' -o matrices/vmtxd/materialMapSouth.hdr -o matrices/vmtxd/hdrIllum/south${idx}.hdr > matrices/vmtxd/hdrLum/south${idx}.hdr
    # done
    # aka `matrices/vmtxd/hdrLum/south${idx}.hdr`
    @staticmethod
    def luminance_images_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("luminance_images",)
            ),
            ShellCommand(
                formatted_command=(
                    "for idx in {{000..144}}\n"
                    "do\n"
                    "pcomb -h -e 'ro=ri(1)*ri(2);go=gi(1)*gi(2);bo=bi(1)*bi(2)'"
                    " -o {material_map_image} -o {images_view_matrix}"
                    " > {output}\n"
                    "done"
                ),
                label_to_input={
                    DataFlowLabel("material_map_image"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,  # Note that this input may probably be streamed but is read 145 times because of the loop, so using files is more appropriate!
                    ),
                    DataFlowLabel("images_view_matrix"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%03d", wildcard="${idx}"),
                    ),
                },
                output=Output(
                    data_format=DataFormat.HDR,
                    series=Series(format="%03d", wildcard="${idx}"),
                ),
            ),
            {
                DataFlowLabel(
                    "material_map_image"
                ): DirectSolarPartOfThreePhaseMethodSimulation.material_map_image_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                ),
                DataFlowLabel(
                    "images_view_matrix"
                ): DirectSolarPartOfThreePhaseMethodSimulation.images_view_matrix_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                ),
            },
        )

    ##D matrix
    # rfluxmtx -v -ff -ab 0 -ad 1000 -lw 0.001 -c 1000 -n 16 objects/GlazingVmtx.rad skyDomes/skyglow.rad -i octrees/room3phDirect.oct > matrices/dmtxd/daylight.dmx
    # aka `matrices/dmtxd/daylight.dmx`
    @staticmethod
    def daylight_transmission_matrix_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        outside_scene_path: Path,
        room_scene_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("daylight_transmission_matrix",)
            ),
            ShellCommand(
                formatted_command="rfluxmtx -v -ff -ab 0 -ad 1000 -lw 0.001 -c 1000 -n 16 {geometry} {scene} -i {octree}",
                label_to_input={
                    DataFlowLabel("geometry"): Input(data_format=DataFormat.RAD),
                    DataFlowLabel("scene"): Input(
                        data_format=DataFormat.RAD, repetition_count=2
                    ),
                    DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                },
                output=Output(data_format=DataFormat.DMX),
            ),
            {
                DataFlowLabel("geometry"): Root.glazing_aperture_geometry_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                ),
                DataFlowLabel("scene"): Root.outside_scene_graph(
                    outside_scene_path=outside_scene_path
                ),
                DataFlowLabel("octree"): Root.Black.black_room_octree_graph(
                    room_scene_path=room_scene_path
                ),
            },
        )

    ##A direct-only sky-matrix
    # gendaymtx -m 1 -d assets/NYC.wea > skyVectors/NYCd.smx
    # aka `skyVectors/NYCd.smx`
    @staticmethod
    def annual_sky_matrix_graph(
        *, typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("annual_sky_matrix",)
            ),
            ShellCommand(
                formatted_command="gendaymtx -m 1 -d {weather}",
                label_to_input={
                    DataFlowLabel("weather"): Input(data_format=DataFormat.WEA)
                },
                output=Output(data_format=DataFormat.SMX),
            ),
            {
                DataFlowLabel(
                    "weather"
                ): Root.typical_meteorological_year_weather_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                )
            },
        )

    ##RESULTS
    ###Illuminance

    # dctimestep matrices/vmtxd/v.mtx matrices/tmtx/blinds.xml matrices/dmtxd/daylight.dmx skyVectors/NYCd.smx
    #  | rmtxop -fa -t -c 47.4 119.9 11.6 -
    #  > results/5ph/3phdir/3phAnnual.ill
    # aka `results/5ph/3phdir/3phAnnual.ill`
    @staticmethod
    def annual_illuminance_graph(
        *,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        outside_scene_path: Path,
        points_path: Path,
        room_scene_path: Path,
        typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH
                + ("annual_illuminance",)
            ),
            ShellCommand(
                formatted_command=(
                    "dctimestep {view_matrix} {bsdf} {transmission_matrix} {sky_matrix}"
                    " | rmtxop -fa -t -c 47.4 119.9 11.6 -"
                ),
                label_to_input={
                    DataFlowLabel("view_matrix"): Input(data_format=DataFormat.MTX),
                    DataFlowLabel("bsdf"): Input(data_format=DataFormat.XML),
                    DataFlowLabel("transmission_matrix"): Input(
                        data_format=DataFormat.DMX
                    ),
                    DataFlowLabel("sky_matrix"): Input(data_format=DataFormat.SMX),
                },
                output=Output(data_format=DataFormat.ILL),
            ),
            {
                DataFlowLabel(
                    "view_matrix"
                ): DirectSolarPartOfThreePhaseMethodSimulation.illuminance_view_matrix_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    points_path=points_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                ),
                DataFlowLabel("bsdf"): Root.bsdf_xml_graph(bsdf_xml_path=bsdf_xml_path),
                DataFlowLabel(
                    "transmission_matrix"
                ): DirectSolarPartOfThreePhaseMethodSimulation.daylight_transmission_matrix_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    outside_scene_path=outside_scene_path,
                    room_scene_path=room_scene_path,
                ),
                DataFlowLabel(
                    "sky_matrix"
                ): DirectSolarPartOfThreePhaseMethodSimulation.annual_sky_matrix_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                ),
            },
        )

    ###Images
    # dctimestep -o results/5ph/3phdir/hdr/south%04d.hdr matrices/vmtxd/hdrLum/south%03d.hdr matrices/tmtx/blinds.xml matrices/dmtxd/daylight.dmx skyVectors/NYCd.smx
    # aka `results/5ph/3phdir/hdr/south%04d.hdr`
    @staticmethod
    def images_graph(
        *,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        outside_scene_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        view_path: Path,
        typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSolarPartOfThreePhaseMethodSimulation.GROUP_PATH + ("images",)
            ),
            ShellCommand(
                formatted_command="dctimestep -o {output} {luminance_images} {bsdf} {transmission_matrix} {sky_matrix}",
                label_to_input={
                    DataFlowLabel("luminance_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%03d"),
                    ),
                    # TODO Check whether the following inputs may be streamed!
                    DataFlowLabel("bsdf"): Input(
                        data_format=DataFormat.XML,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("transmission_matrix"): Input(
                        data_format=DataFormat.DMX,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("sky_matrix"): Input(
                        data_format=DataFormat.SMX,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                },
                output=Output(data_format=DataFormat.HDR, series=Series(format="%04d")),
            ),
            {
                DataFlowLabel(
                    "luminance_images"
                ): DirectSolarPartOfThreePhaseMethodSimulation.luminance_images_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                ),
                DataFlowLabel("bsdf"): Root.bsdf_xml_graph(bsdf_xml_path=bsdf_xml_path),
                DataFlowLabel(
                    "transmission_matrix"
                ): ThreePhaseMethod.Simulation.daylight_transmission_matrix_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    outside_scene_path=outside_scene_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_material_path=room_material_path,
                    room_scene_path=room_scene_path,
                ),
                DataFlowLabel(
                    "sky_matrix"
                ): DirectSolarPartOfThreePhaseMethodSimulation.annual_sky_matrix_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                ),
            },
        )


class DirectSunCoefficients:

    GROUP_PATH = ("direct_sun_coefficients",)
    # PART3 : DIRECT SUN COEFFICIENTS SIMULATION

    ##Create sun primitive definition for solar calculations.
    # echo "void light solar 0 0 3 1e6 1e6 1e6"
    #  > skies/suns.rad
    ##Create solar discs and corresponding modifiers for 5185 suns corresponding to a
    ##Reinhart MF:6 subdivision.
    ##Windows users, who are unlikely to be able to run a MF:6 simulation, should set
    ##"cnt 5185" to "cnt 145" and "MF:6" to "MF:1".
    # cnt 5185 | rcalc -e MF:6 -f reinsrc.cal -e Rbin=recno
    #  -o 'solar source sun 0 0 4 ${Dx} ${Dy} ${Dz} 0.533'
    #  >> skies/suns.rad
    # aka `skies/suns.rad`
    @staticmethod
    def suns_graph() -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(DirectSunCoefficients.GROUP_PATH + ("suns",)),
            ShellCommand(
                formatted_command=(
                    'echo "void light solar 0 0 3 1e6 1e6 1e6"'
                    " && (cnt 5185 | rcalc -e MF:6 -f {calculation} -e Rbin=recno"
                    " -o 'solar source sun 0 0 4 ${{Dx}} ${{Dy}} ${{Dz}} 0.533')"
                ),
                label_to_input={
                    DataFlowLabel("calculation"): Input(data_format=DataFormat.CAL)
                },
                output=Output(data_format=DataFormat.RAD),
            ),
            {
                DataFlowLabel(
                    "calculation"
                ): Root.reinhart_sky_directions_calculation_graph()
            },
        )

    ##Create an octree black octree, shading device with proxy BSDFs and solar discs.
    # oconv -f materialBlack.rad roomBlack.rad skies/suns.rad blinds/blindsWithProxy.rad > octrees/sunCoefficients.oct
    # aka `octrees/sunCoefficients.oct`
    @staticmethod
    def black_shading_device_octree_graph(
        *,
        materialless_glazing_aperture_geometry_path: Path,
        tensor_tree_path: Path,
        room_scene_path: Path,
        glazing_aperture_geometry_material_name: str
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSunCoefficients.GROUP_PATH + ("black_shading_device_octree",)
            ),
            ShellCommand(
                formatted_command="oconv -f {material} {scene} {suns} {bsdf}",
                label_to_input={
                    DataFlowLabel("material"): Input(data_format=DataFormat.MAT),
                    DataFlowLabel("scene"): Input(data_format=DataFormat.RAD),
                    DataFlowLabel("suns"): Input(data_format=DataFormat.RAD),
                    DataFlowLabel("bsdf"): Input(data_format=DataFormat.RAD),
                },
                output=Output(data_format=DataFormat.OCT),
            ),
            {
                DataFlowLabel("material"): Root.Black.black_material_graph(),
                DataFlowLabel("scene"): Root.Black.black_room_scene_graph(
                    room_scene_path=room_scene_path
                ),
                DataFlowLabel("suns"): DirectSunCoefficients.suns_graph(),
                DataFlowLabel("bsdf"): Root.bsdf_with_proxy_rad_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    tensor_tree_path=tensor_tree_path,
                ),
            },
        )

    ##Calculate illuminance sun coefficients for images for the view file south.vf.
    # vwrays -vf views/south.vf -x 400 -y 400 -pj 0.7 -ff
    #   | rcontrib -w- -i -ab 1 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -ffc -n 16
    #   `vwrays -vf views/south.vf -x 400 -y 400 -d`
    #   -o matrices/cds/hdrIllSpace/southM6%04d.hdr -e MF:6 -f reinhart.cal -b rbin -bn Nrbins -m solar octrees/sunCoefficients.oct
    # aka `matrices/cds/hdrIllSpace/southM6%04d.hdr`
    @staticmethod
    def illuminance_images_graph(
        *,
        tensor_tree_path: Path,
        view_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        glazing_aperture_geometry_material_name: str
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSunCoefficients.GROUP_PATH + ("illuminance_images",)
            ),
            ShellCommand(
                formatted_command=(
                    "vwrays -vf {view} -x 400 -y 400 -pj 0.7 -ff"
                    " | rcontrib -i -ab 1 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -ffc -n 16"
                    " `vwrays -vf {view} -x 400 -y 400 -d`"
                    " -o {output} -e MF:6 -f {calculation} -b rbin -bn Nrbins -m solar {octree}"
                ),
                label_to_input={
                    DataFlowLabel("view"): Input(data_format=DataFormat.VF),
                    DataFlowLabel("calculation"): Input(data_format=DataFormat.CAL),
                    DataFlowLabel("octree"): Input(
                        data_format=DataFormat.OCT,
                        transmission_kind=TransmissionKind.FILE,
                    ),  # TODO The oct-file is read once per output. It can probably be a stream though.
                },
                output=Output(data_format=DataFormat.HDR, series=Series(format="%04d")),
            ),
            {
                DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                DataFlowLabel(
                    "calculation"
                ): Root.reinhart_high_density_sky_patches_calculation_graph(),
                DataFlowLabel(
                    "octree"
                ): DirectSunCoefficients.black_shading_device_octree_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    tensor_tree_path=tensor_tree_path,
                    room_scene_path=room_scene_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
            },
        )

    ##Calculate illuminance sun coefficients for illuminance calculations.
    # rcontrib -I+ -ab 1 -y 100 -n 16 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -faf -e MF:6 -f reinhart.cal -b rbin -bn Nrbins -m solar octrees/sunCoefficients.oct < points.txt > matrices/cds/cds.mtx
    # aka `matrices/cds/cds.mtx`
    @staticmethod
    def cds_mtx_graph(
        *,
        materialless_glazing_aperture_geometry_path: Path,
        points_path: Path,
        tensor_tree_path: Path,
        room_scene_path: Path,
        glazing_aperture_geometry_material_name: str
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(DirectSunCoefficients.GROUP_PATH + ("cds_mtx",)),
            ShellCommand(
                formatted_command="rcontrib -I+ -ab 1 -y 100 -n 16 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -faf -e MF:6 -f {calculation} -b rbin -bn Nrbins -m solar {octree} < {points}",
                label_to_input={
                    DataFlowLabel("calculation"): Input(data_format=DataFormat.CAL),
                    # We transmit the octree as file because when running the whole computation for an unknown reason we got the following error:
                    # The command 'set -euo pipefail && (tee /tmp/tmpj6lehvwk/direct_sun_coefficients/pipe_from_black_shading_device_octree_to_direct_sun_coefficients_cds_mtx.oct)'
                    # exited with return code '141' and error message 'b'tee: /tmp/tmpj6lehvwk/direct_sun_coefficients/pipe_from_black_shading_device_octree_to_direct_sun_coefficients_cds_mtx.oct: I/O error\n''
                    DataFlowLabel("octree"): Input(
                        data_format=DataFormat.OCT,
                        # repetition_count=2,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("points"): Input(data_format=DataFormat.PTS),
                },
                output=Output(data_format=DataFormat.MTX),
            ),
            {
                DataFlowLabel(
                    "calculation"
                ): Root.reinhart_high_density_sky_patches_calculation_graph(),
                DataFlowLabel(
                    "octree"
                ): DirectSunCoefficients.black_shading_device_octree_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    tensor_tree_path=tensor_tree_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    room_scene_path=room_scene_path,
                ),
                DataFlowLabel("points"): Root.points_graph(points_path=points_path),
            },
        )

    # aka `objects/GlazingVmtxBlack.rad`
    # data_flow_graph.add_node(
    #     ComputationLabel(
    #         GROUP_PATH + ("black_glazing_aperture_geometry",)
    #     ),
    #     ShellCommand(
    #         formatted_command=(
    #             "echo -e 'void plastic {material_name} 0 0 5 0 0 0 0 0\n'".format(
    #                 material_name=glazing_aperture_geometry_material_name
    #             )
    #             + " && cat {materialless_geometry}"
    #         ),
    #         label_to_input={
    #             DataFlowLabel("materialless_geometry"): Input(data_format=DataFormat.RAD)
    #         },
    #         output=Output(data_format=DataFormat.RAD),
    #     ),
    #     {
    #         DataFlowLabel("materialless_geometry"): ComputationLabel(
    #             ("materialless_glazing_aperture_geometry",)
    #         )
    #     },
    # )

    ##Create a material map
    # oconv -f room.mat room.rad objects/GlazingVmtxBlack.rad > octrees/materialMap5Phase.oct
    # aka `octrees/materialMap5Phase.oct`
    # data_flow_graph.add_node(
    #     ComputationLabel(GROUP_PATH + ("material_map",)),
    #     _frozen_octree_from_material_and_scene_and_geometry(),
    #     {
    #         DataFlowLabel("material"): ComputationLabel(("room_material",)),
    #         DataFlowLabel("scene"): ComputationLabel(("room_scene",)),
    #         DataFlowLabel("geometry"): ComputationLabel(
    #             GROUP_PATH
    #             + ("black_glazing_aperture_geometry",)
    #         ),
    #     },
    # )

    # data_flow_graph.add_node(
    #     ComputationLabel(GROUP_PATH + ("save_material_map",)),
    #     SaveToFile(
    #         label_to_input_with_path={
    #             DataFlowLabel("data"): (
    #                 Input(data_format=DataFormat.OCT),
    #                 temporary_file_path(
    #                     ComputationLabel(
    #                         GROUP_PATH + ("material_map",)
    #                     )
    #                 ),
    #             )
    #         }
    #     ),
    #     {
    #         DataFlowLabel("data"): ComputationLabel(
    #             GROUP_PATH + ("material_map",)
    #         )
    #     },
    # )

    # rpict -x 400 -y 400 -ps 1 -av 0.31831 0.31831 0.31831 -ab 0 -vf views/south.vf octrees/materialMap3PhaseDirect.oct > matrices/cds/materialMapSouth.hdr
    # aka `matrices/cds/materialMapSouth.hdr`
    @staticmethod
    def material_map_images_view_matrix_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSunCoefficients.GROUP_PATH + ("material_map_images_view_matrix",)
            ),
            ShellCommand(
                formatted_command="rpict -x 400 -y 400 -ps 1 -av 0.31831 0.31831 0.31831 -ab 0 -vf {view} {octree}",
                label_to_input={
                    DataFlowLabel("view"): Input(data_format=DataFormat.VF),
                    DataFlowLabel("octree"): Input(data_format=DataFormat.OCT),
                },
                output=Output(data_format=DataFormat.HDR),
            ),
            {
                DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                DataFlowLabel("octree"): Root.Black.opaque_glazing_octree_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                ),
            },
        )

    ##Calculate luminance sun coefficients for images for the view file south.vf.
    # vwrays -vf views/south.vf -x 400 -y 400 -pj 0.7 -ff
    #   | rcontrib -w- -ab 1 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -ffc -n 16
    #   `vwrays -vf views/south.vf -x 400 -y 400 -d`
    #   -o matrices/cds/hdrLumFacade/southM6%04d.hdr -e MF:6 -f reinhart.cal -b rbin -bn Nrbins -m solar octrees/sunCoefficients.oct
    # aka `matrices/cds/hdrLumFacade/southM6%04d.hdr`
    @staticmethod
    def luminance_images_graph(
        *,
        tensor_tree_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        view_path: Path,
        glazing_aperture_geometry_material_name: str
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(DirectSunCoefficients.GROUP_PATH + ("luminance_images",)),
            ShellCommand(
                formatted_command=(
                    "vwrays -vf {view} -x 400 -y 400 -pj 0.7 -ff"
                    " | rcontrib -ab 1 -ad 256 -lw 1.0e-3 -dc 1 -dt 0 -dj 0 -ffc -n 16"
                    " `vwrays -vf {view} -x 400 -y 400 -d`"
                    " -o {output} -e MF:6 -f {calculation} -b rbin -bn Nrbins -m solar {octree}"
                ),
                label_to_input={
                    DataFlowLabel("view"): Input(data_format=DataFormat.VF),
                    DataFlowLabel("calculation"): Input(data_format=DataFormat.CAL),
                    DataFlowLabel("octree"): Input(
                        data_format=DataFormat.OCT,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                },
                output=Output(data_format=DataFormat.HDR, series=Series(format="%04d")),
            ),
            {
                DataFlowLabel("view"): Root.view_graph(view_path=view_path),
                DataFlowLabel(
                    "calculation"
                ): Root.reinhart_high_density_sky_patches_calculation_graph(),
                DataFlowLabel(
                    "octree"
                ): DirectSunCoefficients.black_shading_device_octree_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    tensor_tree_path=tensor_tree_path,
                    room_scene_path=room_scene_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
            },
        )

    ##Calculate luminance sun coefficients for images for the view file views/south.vf.
    # for idx in {0000..5185}
    # do
    #     pcomb -h -e 'ro=ri(1)*ri(2)+ri(3);go=gi(1)*gi(2)+gi(3);bo=bi(1)*bi(2)+bi(3)' -o matrices/cds/materialMapSouth.hdr -o matrices/cds/hdrIllSpace/southM6${idx}.hdr -o matrices/cds/hdrLumFacade/southM6${idx}.hdr > matrices/cds/hdr/south${idx}.hdr
    # done
    # aka `matrices/cds/hdr/south${idx}.hdr`
    @staticmethod
    def luminance_images_2_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSunCoefficients.GROUP_PATH + ("luminance_images_2",)
            ),
            ShellCommand(
                formatted_command=(
                    "for idx in {{0000..5185}}\n"
                    "do\n"
                    "pcomb -h -e 'ro=ri(1)*ri(2)+ri(3);go=gi(1)*gi(2)+gi(3);bo=bi(1)*bi(2)+bi(3)' -o {material_map} -o {illuminance_images} -o {luminance_images} > {output}\n"
                    "done"
                ),
                label_to_input={
                    DataFlowLabel("material_map"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,  # Note that this input may probably be streamed but is read 5186 times because of the loop, so using files is more appropriate!
                    ),
                    DataFlowLabel("illuminance_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d", wildcard="${idx}"),
                    ),
                    DataFlowLabel("luminance_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d", wildcard="${idx}"),
                    ),
                },
                output=Output(
                    data_format=DataFormat.HDR,
                    series=Series(format="%04d", wildcard="${idx}"),
                ),
            ),
            {
                DataFlowLabel(
                    "material_map"
                ): DirectSunCoefficients.material_map_images_view_matrix_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                ),
                DataFlowLabel(
                    "illuminance_images"
                ): DirectSunCoefficients.illuminance_images_graph(
                    tensor_tree_path=tensor_tree_path,
                    view_path=view_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
                DataFlowLabel(
                    "luminance_images"
                ): DirectSunCoefficients.luminance_images_graph(
                    tensor_tree_path=tensor_tree_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
            },
        )

    ##Generate a sun-matrix for the sun coefficients calculation.
    # gendaymtx -5 0.533 -d -m 6 assets/NYC.wea > skyVectors/NYCsun.smx
    # aka `skyVectors/NYCsun.smx`
    @staticmethod
    def sun_matrix_graph(
        *, typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(DirectSunCoefficients.GROUP_PATH + ("sun_matrix",)),
            ShellCommand(
                formatted_command="gendaymtx -5 0.533 -d -m 6 {weather}",
                label_to_input={
                    DataFlowLabel("weather"): Input(data_format=DataFormat.WEA)
                },
                output=Output(data_format=DataFormat.SMX),
            ),
            {
                DataFlowLabel(
                    "weather"
                ): Root.typical_meteorological_year_weather_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                )
            },
        )

    ###Images
    # dctimestep -o results/5ph/cds/hdr/south%04d.hdr matrices/cds/hdr/south%04d.hdr skyVectors/NYCsun.smx
    # aka `results/5ph/cds/hdr/south%04d.hdr`
    @staticmethod
    def images_graph(
        *,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        typical_meteorological_year_weather_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(DirectSunCoefficients.GROUP_PATH + ("images",)),
            ShellCommand(
                formatted_command="dctimestep -o {output} {luminance_images} {sun_matrix}",
                label_to_input={
                    DataFlowLabel("luminance_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d"),
                    ),
                    DataFlowLabel("sun_matrix"): Input(data_format=DataFormat.SMX),
                },
                output=Output(data_format=DataFormat.HDR, series=Series(format="%04d")),
            ),
            {
                DataFlowLabel(
                    "luminance_images"
                ): DirectSunCoefficients.luminance_images_2_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_scene_path=room_scene_path,
                    tensor_tree_path=tensor_tree_path,
                    view_path=view_path,
                ),
                DataFlowLabel("sun_matrix"): DirectSunCoefficients.sun_matrix_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                ),
            },
        )

    ###Illuminance
    # dctimestep matrices/cds/cds.mtx skyVectors/NYCsun.smx
    #  | rmtxop -fa -t -c 47.4 119.9 11.6 - > results/5ph/cds/cds.ill
    # aka `results/5ph/cds/cds.ill`
    @staticmethod
    def annual_illuminance_graph(
        *,
        materialless_glazing_aperture_geometry_path: Path,
        typical_meteorological_year_weather_path: Path,
        points_path: Path,
        tensor_tree_path: Path,
        room_scene_path: Path,
        glazing_aperture_geometry_material_name: str
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(
                DirectSunCoefficients.GROUP_PATH + ("annual_illuminance",)
            ),
            ShellCommand(
                formatted_command=(
                    "dctimestep {cds_mtx} {sun_matrix} | rmtxop -fa -t -c 47.4 119.9 11.6 -"
                ),
                label_to_input={
                    DataFlowLabel("cds_mtx"): Input(
                        data_format=DataFormat.MTX,
                        transmission_kind=TransmissionKind.FILE,
                    ),
                    DataFlowLabel("sun_matrix"): Input(data_format=DataFormat.SMX),
                },
                output=Output(data_format=DataFormat.ILL),
            ),
            {
                DataFlowLabel("cds_mtx"): DirectSunCoefficients.cds_mtx_graph(
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    points_path=points_path,
                    tensor_tree_path=tensor_tree_path,
                    room_scene_path=room_scene_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
                DataFlowLabel("sun_matrix"): DirectSunCoefficients.sun_matrix_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path
                ),
            },
        )


class FivePhaseMethod:

    GROUP_PATH = ("five_phase_method",)
    ##Combine Image results using the 5 Phase Equation.
    # for idx in {0001..0040}
    # do
    #     pcomb -h -e 'ro=ri(1)-ri(2)+ri(3);go=gi(1)-gi(2)+gi(3);bo=bi(1)-bi(2)+bi(3)' -o results/5ph/3ph/hdr/south${idx}.hdr -o results/5ph/3phdir/hdr/south${idx}.hdr -o -o results/5ph/cds/hdr/south${idx}.hdr > results/5ph/hdr/south${idx}.hdr
    # done
    # aka `results/5ph/hdr/south${idx}.hdr`
    @staticmethod
    def images_graph(
        *,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        outside_scene_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        typical_meteorological_year_weather_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(FivePhaseMethod.GROUP_PATH + ("images",)),
            ShellCommand(
                formatted_command=(
                    "for idx in {{0001..0040}}\n"
                    "do\n"
                    "pcomb -h"
                    " -e 'ro=ri(1)-ri(2)+ri(3);go=gi(1)-gi(2)+gi(3);bo=bi(1)-bi(2)+bi(3)'"
                    " -o {three_phase_method_images}"
                    " -o {direct_solar_part_images}"
                    " -o"
                    " -o {direct_sun_coefficients_images}"
                    " > {output}\n"
                    "done"
                ),
                label_to_input={
                    DataFlowLabel("three_phase_method_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d", wildcard="${idx}"),
                    ),
                    DataFlowLabel("direct_solar_part_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d", wildcard="${idx}"),
                    ),
                    DataFlowLabel("direct_sun_coefficients_images"): Input(
                        data_format=DataFormat.HDR,
                        transmission_kind=TransmissionKind.FILE,
                        series=Series(format="%04d", wildcard="${idx}"),
                    ),
                },
                output=Output(
                    data_format=DataFormat.HDR,
                    series=Series(format="%04d", wildcard="${idx}"),
                ),
            ),
            {
                DataFlowLabel(
                    "three_phase_method_images"
                ): ThreePhaseMethod.Simulation.images_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    outside_scene_path=outside_scene_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    room_material_path=room_material_path,
                    room_scene_path=room_scene_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                    view_path=view_path,
                ),
                DataFlowLabel(
                    "direct_solar_part_images"
                ): DirectSolarPartOfThreePhaseMethodSimulation.images_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    outside_scene_path=outside_scene_path,
                    room_material_path=room_material_path,
                    room_scene_path=room_scene_path,
                    view_path=view_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                ),
                DataFlowLabel(
                    "direct_sun_coefficients_images"
                ): DirectSunCoefficients.images_graph(
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                    room_scene_path=room_scene_path,
                    tensor_tree_path=tensor_tree_path,
                    view_path=view_path,
                ),
            },
        )

    @staticmethod
    def save_images_graph(
        *,
        images_output_path: Path,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        outside_scene_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        typical_meteorological_year_weather_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        """Returns the data-flow graph that when executed computes the images
        with the five-phase method for the given inputs."""
        return DataFlowGraph.construct(
            ComputationLabel(FivePhaseMethod.GROUP_PATH + ("save_images",)),
            SaveToFile(
                label_to_input_with_path={
                    DataFlowLabel("data"): (
                        Input(
                            data_format=DataFormat.HDR,
                            transmission_kind=TransmissionKind.FILE,
                            series=Series(format="%04d", wildcard="*"),
                        ),
                        images_output_path,
                    )
                }
            ),
            {
                DataFlowLabel("data"): FivePhaseMethod.images_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    outside_scene_path=outside_scene_path,
                    room_material_path=room_material_path,
                    room_scene_path=room_scene_path,
                    tensor_tree_path=tensor_tree_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                    view_path=view_path,
                )
            },
        )

    # Combine illuminance results using the 5 Phase equation.
    # rmtxop results/5ph/3ph/3phAnnual.ill + -s -1 results/5ph/3phdir/3phAnnual.ill + results/5ph/cds/cds.ill > results/5ph/5Phannual.ill
    # aka `results/5ph/5Phannual.ill`
    # The `sed` commands remove the preamble, unnecessary white-space, and
    # make the matrix comma separated instead of tab separated.
    @staticmethod
    def annual_illuminance_graph(
        *,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        materialless_glazing_aperture_geometry_path: Path,
        outside_scene_path: Path,
        points_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        return DataFlowGraph.construct(
            ComputationLabel(FivePhaseMethod.GROUP_PATH + ("annual_illuminance",)),
            ShellCommand(
                formatted_command=(  # pylint: disable=anomalous-backslash-in-string
                    "rmtxop {three_phase_method_annual_illuminance} + -s -1 {direct_solar_part_annual_illuminance} + {direct_sun_coefficients_annual_illuminance}"
                    " | sed -e '1,/^$/d'"  # Remove all lines until and including the first blank line
                    " | sed -e 's/^[ \t]\+\|[ \t]\+$//g'"  # Remove spaces and tabs at beginning and end of lines
                    " | sed -e 's/[ \t]\+/,/g'"  # Replace spaces and tabs by commas
                ),
                label_to_input={
                    DataFlowLabel("three_phase_method_annual_illuminance"): Input(
                        data_format=DataFormat.ILL
                    ),
                    DataFlowLabel("direct_solar_part_annual_illuminance"): Input(
                        data_format=DataFormat.ILL
                    ),
                    DataFlowLabel("direct_sun_coefficients_annual_illuminance"): Input(
                        data_format=DataFormat.ILL
                    ),
                },
                output=Output(data_format=DataFormat.ILL),
            ),
            {
                DataFlowLabel(
                    "three_phase_method_annual_illuminance"
                ): ThreePhaseMethod.Simulation.annual_illuminance_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    points_path=points_path,
                    outside_scene_path=outside_scene_path,
                    room_material_path=room_material_path,
                    room_scene_path=room_scene_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                ),
                DataFlowLabel(
                    "direct_solar_part_annual_illuminance"
                ): DirectSolarPartOfThreePhaseMethodSimulation.annual_illuminance_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    points_path=points_path,
                    outside_scene_path=outside_scene_path,
                    room_scene_path=room_scene_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                ),
                DataFlowLabel(
                    "direct_sun_coefficients_annual_illuminance"
                ): DirectSunCoefficients.annual_illuminance_graph(
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    points_path=points_path,
                    tensor_tree_path=tensor_tree_path,
                    room_scene_path=room_scene_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                ),
            },
        )

    @staticmethod
    def save_annual_illuminance_graph(
        *,
        annual_illuminance_output_path: Path,
        bsdf_xml_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        glazing_aperture_geometry_material_name: str,
        outside_scene_path: Path,
        points_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        tensor_tree_path: Path,
        typical_meteorological_year_weather_path: Path
    ) -> DataFlowGraph:
        """Returns the data-flow graph that when executed computes the annual
        illuminance with the five-phase method for the given inputs."""
        return DataFlowGraph.construct(
            ComputationLabel(FivePhaseMethod.GROUP_PATH + ("save_annual_illuminance",)),
            SaveToFile(
                label_to_input_with_path={
                    DataFlowLabel("data"): (
                        Input(data_format=DataFormat.ILL),
                        annual_illuminance_output_path,
                    )
                }
            ),
            {
                DataFlowLabel("data"): FivePhaseMethod.annual_illuminance_graph(
                    bsdf_xml_path=bsdf_xml_path,
                    glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
                    materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
                    outside_scene_path=outside_scene_path,
                    points_path=points_path,
                    room_scene_path=room_scene_path,
                    typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
                    room_material_path=room_material_path,
                    tensor_tree_path=tensor_tree_path,
                )
            },
        )

    @staticmethod
    def graph(
        *,
        annual_illuminance_output_path: Path,
        bsdf_xml_path: Path,
        glazing_aperture_geometry_material_name: str,
        images_output_path: Path,
        materialless_glazing_aperture_geometry_path: Path,
        points_path: Path,
        room_material_path: Path,
        room_scene_path: Path,
        outside_scene_path: Path,
        tensor_tree_path: Path,
        typical_meteorological_year_weather_path: Path,
        view_path: Path
    ) -> DataFlowGraph:
        """Returns the data-flow graph that when executed computes the annual
        illuminance and the images with the five-phase method for the given
        inputs."""
        save_annual_illuminance_graph = FivePhaseMethod.save_annual_illuminance_graph(
            annual_illuminance_output_path=annual_illuminance_output_path,
            bsdf_xml_path=bsdf_xml_path,
            glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
            materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
            outside_scene_path=outside_scene_path,
            points_path=points_path,
            room_scene_path=room_scene_path,
            typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
            room_material_path=room_material_path,
            tensor_tree_path=tensor_tree_path,
        )
        save_images_graph = FivePhaseMethod.save_images_graph(
            images_output_path=images_output_path,
            bsdf_xml_path=bsdf_xml_path,
            glazing_aperture_geometry_material_name=glazing_aperture_geometry_material_name,
            materialless_glazing_aperture_geometry_path=materialless_glazing_aperture_geometry_path,
            outside_scene_path=outside_scene_path,
            room_material_path=room_material_path,
            room_scene_path=room_scene_path,
            tensor_tree_path=tensor_tree_path,
            typical_meteorological_year_weather_path=typical_meteorological_year_weather_path,
            view_path=view_path,
        )
        return LabeledDirectedAcyclicGraph.merge(
            [save_annual_illuminance_graph, save_images_graph]
        )


if __name__ == "__main__":
    compute(
        annual_illuminance_output_path=Path(
            "tmp", *FivePhaseMethod.GROUP_PATH, "annual_illuminance"
        ),
        bsdf_xml_path=Path(
            "tests/test_five_phase_method/test_integration/matrices/tmtx/blinds.xml"
        ),
        glazing_aperture_geometry_material_name="Glazing",
        images_output_path=Path("tmp", *FivePhaseMethod.GROUP_PATH, "images"),
        materialless_glazing_aperture_geometry_path=Path(
            "tests/test_five_phase_method/test_integration/objects/GlazingVmtx.rad"
        ),
        outside_scene_path=Path(
            "tests/test_five_phase_method/test_integration/skyDomes/skyglow.rad"
        ),
        points_path=Path("tests/test_five_phase_method/test_integration/points.pts"),
        room_material_path=Path(
            "tests/test_five_phase_method/test_integration/room.mat"
        ),
        room_scene_path=Path("tests/test_five_phase_method/test_integration/room.rad"),
        tensor_tree_path=Path(
            "tests/test_five_phase_method/test_integration/blinds/blinds_t4c.xml"
        ),
        typical_meteorological_year_weather_path=Path(
            "tests/test_five_phase_method/test_integration/assets/NYC.wea"
        ),
        view_path=Path("tests/test_five_phase_method/test_integration/views/south.vf"),
    )
    # pass
