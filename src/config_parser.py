from configparser import ConfigParser
from numpy import ndarray


class Config:
    pass


def parse(config_path, options):
    parser = ConfigParser()
    parser.read(config_path)
    config = Config()
    _add_simulation_variables(parser, config)  # adds `numConWin`
    _add_geometry_variables(parser, config, options)
    _add_lighting_schedule_variables(parser, config, options)
    _add_daylighting_variables(parser, config, options, config.numConWin)
    _add_glare_variables(parser, config, options, config.numConWin)
    _add_thermal_variables(parser, config, options, config.numConWin)
    _add_cut_off_variables(parser, config, options)
    _add_matrix_control_variables(parser, config, options)
    _add_shading_control_variables(parser, config, options)
    _add_klems_variables(parser, config, options, config.numConWin)
    _add_calorim_variables(parser, config, options, config.numConWin)
    _add_ep_variables(parser, config, options, config.numConWin)
    return config


def _add_simulation_variables(parser, config):
    config.meteo_path = parser.get("PATHS", "meteo")
    config.lat = parser.getfloat("VARIABLES", "lat")
    config.lon = parser.getfloat("VARIABLES", "lon")
    config.tzone = parser.getfloat("VARIABLES", "tzone")
    config.iniMonth = parser.getint("VARIABLES", "iniMonth")
    config.iniDay = parser.getint("VARIABLES", "iniDay")
    config.endMonth = parser.getint("VARIABLES", "endMonth")
    config.endDay = parser.getint("VARIABLES", "endDay")
    config.iniDayWeek = parser.getint("VARIABLES", "iniDayWeek")
    config.floor = parser.getint("VARIABLES", "floor")
    config.grndAlb = parser.getfloat("VARIABLES", "grndAlb")
    config.floorArea = parser.getfloat("VARIABLES", "floorArea")
    config.numConWin = parser.getint("VARIABLES", "numConWin")


def _add_geometry_variables(parser, config, options):
    config.surf_path = parser.get("PATHS", "surf")
    config.win_path = parser.get("PATHS", "win")
    config.frame_path = parser.get("PATHS", "frame")
    config.output_path = parser.get("PATHS", "output")
    config.input_path = parser.get("PATHS", "input")
    config.workDir_path = parser.get("PATHS", "workDir")
    if options.shoeBox or options.geo or options.ep:
        config.height = parser.getfloat("VARIABLES", "height")
        config.rotAng = parser.getfloat("VARIABLES", "rotAng")
    if options.shoeBox or options.ep:
        config.length = parser.getfloat("VARIABLES", "length")
        config.width = parser.getfloat("VARIABLES", "width")
        config.thickSouth = parser.getfloat("VARIABLES", "thickSouth")
        config.thickEast = parser.getfloat("VARIABLES", "thickEast")
        config.thickNorth = parser.getfloat("VARIABLES", "thickNorth")
        config.thickWest = parser.getfloat("VARIABLES", "thickWest")
        config.thickCeiling = parser.getfloat("VARIABLES", "thickCeiling")
        config.thickFloor = parser.getfloat("VARIABLES", "thickFloor")
        config.albWall = parser.getfloat("VARIABLES", "albWall")
        config.albCeiling = parser.getfloat("VARIABLES", "albCeiling")
        config.albFloor = parser.getfloat("VARIABLES", "albFloor")
        config.bcSouth = parser.getint("VARIABLES", "bcSouth")
        config.bcEast = parser.getint("VARIABLES", "bcEast")
        config.bcNorth = parser.getint("VARIABLES", "bcNorth")
        config.bcWest = parser.getint("VARIABLES", "bcWest")
        config.bcCeiling = parser.getint("VARIABLES", "bcCeiling")
        config.bcFloor = parser.getint("VARIABLES", "bcFloor")
        config.conSouth = parser.getint("VARIABLES", "conSouth")
        config.conEast = parser.getint("VARIABLES", "conEast")
        config.conNorth = parser.getint("VARIABLES", "conNorth")
        config.conWest = parser.getint("VARIABLES", "conWest")
        config.conCeiling = parser.getint("VARIABLES", "conCeiling")
        config.conFloor = parser.getint("VARIABLES", "conFloor")


def _add_lighting_schedule_variables(parser, config, options):
    if options.lightSch or options.therm or options.ep:
        config.lightSch_path = parser.get("PATHS", "lightSch")
        config.powerLight = parser.getfloat("VARIABLES", "powerLight")
    if options.lightSch:
        config.lightCntrl_path = parser.get("PATHS", "lightCntrl")


def _add_daylighting_variables(parser, config, options, numConWin):
    if (
        options.daylight
        or options.df
        or options.rad == "radBsdf"
        or options.rad == "radTpm"
        or options.therm
        or options.ep
        or options.calorim
        or options.glare
        or options.glareSimpl
        or options.glareMulti
        or options.fpm
    ):
        config.bsdfSys_path = ndarray((numConWin,), dtype=object)
        for i in range(numConWin):
            config.bsdfSys_path[i] = parser.get("PATHS", "bsdfSys_%i" % i)
        if options.fpm:
            config.tensorTree_path = ndarray((numConWin,), dtype=object)
            for i in range(numConWin):
                config.tensorTree_path[i] = parser.get("PATHS", "tensorTree_%i" % i)
    if options.daylight or options.df or options.rad == "radTpm":
        if options.grid:
            config.numPhotosensX = parser.getint("VARIABLES", "numPhotosensX")
            config.numPhotosensY = parser.getint("VARIABLES", "numPhotosensY")
            config.photosensHeight = parser.getfloat("VARIABLES", "photosensHeight")
            config.gridXOffset = parser.getfloat("VARIABLES", "gridXOffset")
            config.gridYOffset = parser.getfloat("VARIABLES", "gridYOffset")
        else:
            config.illuPts_path = parser.get("PATHS", "illuPts")
    if options.da:
        config.numPhotosensX = parser.getint("VARIABLES", "numPhotosensX")
        config.numPhotosensY = parser.getint("VARIABLES", "numPhotosensY")
        config.startHour = parser.getint("VARIABLES", "startHour")
        config.endHour = parser.getint("VARIABLES", "endHour")
        config.minIllum = parser.getint("VARIABLES", "minIllum")
    if options.mask:
        config.obstMask_path = parser.get("PATHS", "obstMask")
    if options.outside:
        config.outside_path = parser.get("PATHS", "outside")


def _add_glare_variables(parser, config, options, numConWin):
    if options.glare or options.glareSimpl or options.glareFull or options.glareMulti:
        config.illuVertPts_path = parser.get("PATHS", "illuVertPts")
    if options.glare or options.glareFull or options.glareMulti:
        config.winRad_path = ndarray((numConWin,), dtype=object)
        for i in range(numConWin):
            config.winRad_path[i] = parser.get("PATHS", "winRad_%i" % i)


def _add_thermal_variables(parser, config, options, numConWin):
    if options.therm or options.ep or options.iso:
        config.matOpaq_path = parser.get("PATHS", "matOpaq")
        config.constOpaq_path = parser.get("PATHS", "constOpaq")
        config.infSch_path = parser.get("PATHS", "infSch")
        config.occSch_path = parser.get("PATHS", "occSch")
        config.equipSch_path = parser.get("PATHS", "equipSch")
        config.heatSpSch_path = parser.get("PATHS", "heatSpSch")
        config.coolSpSch_path = parser.get("PATHS", "coolSpSch")
        config.cop = parser.getfloat("VARIABLES", "cop")
        config.airExch = parser.getfloat("VARIABLES", "airExch")
        config.volume = parser.getfloat("VARIABLES", "volume")
        config.iniTemp = parser.getfloat("VARIABLES", "iniTemp")
        config.powerEquip = parser.getfloat("VARIABLES", "powerEquip")
        config.radFracEquip = parser.getfloat("VARIABLES", "radFracEquip")
        config.powerPeople = parser.getfloat("VARIABLES", "powerPeople")
        config.radFracPeople = parser.getfloat("VARIABLES", "radFracPeople")
        config.radFracLight = parser.getfloat("VARIABLES", "radFracLight")
        config.radFracPeople = parser.getfloat("VARIABLES", "radFracPeople")
        config.calorim_path = ndarray((numConWin,), dtype=object)
        for i in range(numConWin):
            config.calorim_path[i] = parser.get("PATHS", "calorim_%i" % i)
        if options.thermComf:
            config.height = parser.getfloat("VARIABLES", "height")
            config.length = parser.getfloat("VARIABLES", "length")
            config.width = parser.getfloat("VARIABLES", "width")
            config.thickSouth = parser.getfloat("VARIABLES", "thickSouth")
            config.thickEast = parser.getfloat("VARIABLES", "thickEast")
            config.thickNorth = parser.getfloat("VARIABLES", "thickNorth")
            config.thickWest = parser.getfloat("VARIABLES", "thickWest")
            config.thickCeiling = parser.getfloat("VARIABLES", "thickCeiling")
            config.thickFloor = parser.getfloat("VARIABLES", "thickFloor")
    else:
        config.iniTemp = 300


def _add_cut_off_variables(parser, config, options):
    if options.cutOff:
        config.slatWidth = parser.getfloat("VARIABLES", "slatWidth")
        config.slatDist = parser.getfloat("VARIABLES", "slatDist")
        config.minAng = parser.getfloat("VARIABLES", "minAng")
        config.conMinAng = parser.getint("VARIABLES", "conMinAng")
        config.maxAng = parser.getfloat("VARIABLES", "maxAng")
        config.conMaxAng = parser.getint("VARIABLES", "conMaxAng")
        config.stepAng = parser.getfloat("VARIABLES", "stepAng")


def _add_matrix_control_variables(parser, config, options):
    if options.mtxCntrl or options.refeedCntrl or options.optStateCntrl:
        if options.refeedCntrl:
            print("WARNING: control algorithm will refeed according to window 0")
        if options.optStateCntrl:
            print("WARNING: for the optStateCntrl option, windows have to be ordered conIndex = [lowerPartition_win0, lowerPartition_win1,..., upperPartition_win0, upperPartition_win1,...]")
            config.numWinLowerPartition = parser.getint("VARIABLES", "numWinLowerPartition")
            config.numWinUpperPartition = parser.getint("VARIABLES", "numWinUpperPartition")
            config.cntrlOpt = parser.getint("VARIABLES", "cntrlOpt")
        config.cntrlMtx_path = parser.get("PATHS", "cntrlMtx")
        config.conDef = parser.getint("VARIABLES", "conDef")


def _add_shading_control_variables(parser, config, options):
    if options.schCntrl:
        config.shadingSch_path = parser.get("PATHS", "shadingSch")


def _add_klems_variables(parser, config, options, numConWin):
    if options.klems:
        config.matGlz_path = parser.get("PATHS", "matGlz")
        config.constWin_path = parser.get("PATHS", "constWin")
        config.bsdfSys_path = ndarray((numConWin,), dtype=object)
        for i in range(numConWin):
            config.bsdfSys_path[i] = parser.get("PATHS", "bsdfSys_%i" % i)
        config.numBsdfLay = parser.getint("VARIABLES", "numBsdfLay")
        if config.numBsdfLay > 0:
            config.bsdfLay_file = ndarray((config.numBsdfLay,), dtype=object)
            for i in range(config.numBsdfLay):
                config.bsdfLay_file[i] = parser.get("PATHS", "bsdfLay_file_%i" % i)
        config.numPanes = parser.getint("VARIABLES", "numPanes")
        config.absFront_path = ndarray([numConWin, config.numPanes], dtype=object)
        config.absBack_path = ndarray([numConWin, config.numPanes], dtype=object)
        for i in range(numConWin):
            for j in range(config.numPanes):
                config.absFront_path[i, j] = parser.get(
                    "PATHS", "absFront_%i_%i" % (i, j)
                )
                config.absBack_path[i, j] = parser.get(
                    "PATHS", "absBack_%i_%i" % (i, j)
                )


def _add_calorim_variables(parser, config, options, numConWin):
    if options.calorim:
        config.calorim_path = ndarray((numConWin,), dtype=object)
        for i in range(numConWin):
            config.calorim_path[i] = parser.get("PATHS", "calorim_%i" % i)


def _add_ep_variables(parser, config, options, numConWin):
    if options.ep or options.iso or options.calorim:
        config.matGlz_path = parser.get("PATHS", "matGlz")
        config.constWin_path = parser.get("PATHS", "constWin")
        config.matGas_path = parser.get("PATHS", "matGas")
        config.matBSDF_path = parser.get("PATHS", "matBSDF")
        config.numPanes = parser.getint("VARIABLES", "numPanes")
        config.absFront_path = ndarray([numConWin, config.numPanes], dtype=object)
        config.absBack_path = ndarray([numConWin, config.numPanes], dtype=object)
        for i in range(numConWin):
            for j in range(config.numPanes):
                config.absFront_path[i, j] = parser.get(
                    "PATHS", "absFront_%i_%i" % (i, j)
                )
                config.absBack_path[i, j] = parser.get(
                    "PATHS", "absBack_%i_%i" % (i, j)
                )
    if options.ep:
        config.ep_path = parser.get("PATHS", "ep")
        config.frameWidth = parser.getfloat("VARIABLES", "frameWidth")
    if options.iso:
        config.stj_type = parser.getfloat("VARIABLES", "stj_type")
        config.stjTfluidIn = parser.getfloat("VARIABLES", "stjTfluidIn")
        config.stjAref = parser.getfloat("VARIABLES", "stjAref")
        config.stjMassFlow = parser.getfloat("VARIABLES", "stjMassFlow")
        config.stjFluidCp = parser.getfloat("VARIABLES", "stjFluidCp")
