import argparse
import os


def parse():
    return _build_parser().parse_args()


def _build_parser():
    parser = argparse.ArgumentParser()
    _add_common_arguments(parser)
    _add_radiance_subparsers(parser)
    return parser


def _add_common_arguments(parser):
    parser.add_argument("config", help="config file for Fener simulation")
    parser.add_argument(
        "-c",
        action="store",
        type=int,
        default=os.cpu_count() if os.cpu_count() is not None else 1,
        help="number of cores (default: '%(default)s')",
    )
    parser.add_argument(
        "-meteo",
        action="store_true",
        default=False,
        help="generates sky matrix and other weather files",
    )
    parser.add_argument(
        "-geo",
        action="store_true",
        default=False,
        help="generates Radiance geometry and three-phase method matrices",
    )
    parser.add_argument(
        "-shoeBox",
        action="store_true",
        default=False,
        help="generates the file 'surfaces.dat' with a shoe-box geometry. together with the -therm option, it also generates view factors based on geometrical analytical relationships (Mills)",
    )
    parser.add_argument(
        "-grid",
        action="store_true",
        default=False,
        help="generates a grid of illuminance sensors",
    )
    parser.add_argument(
        "-display", action="store_true", default=False, help="displays the geometry"
    )
    parser.add_argument(
        "-daylight",
        action="store_true",
        default=False,
        help="runs a daylight simulation",
    )
    parser.add_argument(
        "-lightSch",
        action="store_true",
        default=False,
        help="calculates a schedule of artificial lighting operation from the daylight simulation to be read in the thermal simulation",
    )
    parser.add_argument(
        "-glare",
        action="store_true",
        default=False,
        help="runs glare simulations based on the Enhanced Simplified Method",
    )
    parser.add_argument(
        "-glareSimpl",
        action="store_true",
        default=False,
        help="runs glare simulations based on the Simplified method",
    )
    parser.add_argument(
        "-glareFull",
        action="store_true",
        default=False,
        help="runs glare simulations based on the full Evalglare method",
    )
    parser.add_argument(
        "-therm", action="store_true", default=False, help="runs thermal simulation"
    )
    parser.add_argument(
        "-thermComf",
        action="store_true",
        default=False,
        help="runs thermal comfort simulation",
    )
    parser.add_argument(
        "-cutOff",
        action="store_true",
        default=False,
        help="includes a cut-off shading control algorithm for blinds",
    )
    parser.add_argument(
        "-mtxCntrl",
        action="store_true",
        default=False,
        help="includes a matrix shading control",
    )
    parser.add_argument(
        "-refeedCntrl",
        action="store_true",
        default=False,
        help="includes a matrix shading control with refeed capability",
    )
    parser.add_argument(
        "-schCntrl",
        action="store_true",
        default=False,
        help="includes a scheduled shading control",
    )
    parser.add_argument(
        "-optStateCntrl",
        action="store_true",
        default=False,
        help="calculates the optimum state for illHor > illHorSp and illVer < illVerSp",
    )
    parser.add_argument(
        "-df",
        action="store_true",
        default=False,
        help="includes a daylight factor calculation",
    )
    parser.add_argument(
        "-da",
        action="store_true",
        default=False,
        help="includes a daylight autonomy and annual sunlight exposure calculation",
    )
    parser.add_argument(
        "-klems",
        action="store_true",
        default=False,
        help="generates a BSDF of a fenestration system from BSDF of the individual layers. For glazing and diffusive layers (shades), the model calculates the BSDF internally",
    )
    parser.add_argument(
        "-mask", action="store_true", default=False, help="generates a obstruction mask"
    )
    parser.add_argument(
        "-outside",
        action="store_true",
        default=False,
        help="generates an outdoor scene",
    )


def _add_radiance_subparsers(parser):
    subparsers = parser.add_subparsers(
        help="run the problem directly with Radiance", dest="rad"
    )
    bsdf_parser = subparsers.add_parser(
        "radBsdf", help="run classic Radiance with BSDF material"
    )
    _add_radiance_bsdf_arguments(bsdf_parser)
    three_phase_method_parser = subparsers.add_parser(
        "radTpm", help="run three-phase method from Radiance"
    )
    _add_radiance_three_phase_method_arguments(three_phase_method_parser)


def _add_radiance_bsdf_arguments(parser):
    parser.add_argument("hour", nargs="+", help="hour of simulation", type=float)
    parser.add_argument(
        "-ab", help="ambient bounces (default: '%(default)s')", type=int, default=4
    )
    parser.add_argument(
        "-ad", help="ambient divisions (default: '%(default)s')", type=int, default=512
    )
    parser.add_argument(
        "-ar", help="ambient resolution (default: '%(default)s')", type=int, default=128
    )
    parser.add_argument(
        "-ass",
        help="ambient super-samples (default: '%(default)s')",
        type=int,
        default=256,
    )
    parser.add_argument(
        "-aa",
        help="ambient accuracy (default: '%(default)s')",
        type=float,
        default=0.15,
    )
    parser.add_argument(
        "-skipDays",
        help="number of days to be skipped from start day (default: '%(default)s')",
        type=int,
        default=0,
    )


def _add_radiance_three_phase_method_arguments(_parser):
    pass
