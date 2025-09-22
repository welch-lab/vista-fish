#!/usr/bin/env jython
# -*- coding: utf-8 -*-

import sys
import os

from ij import IJ
from ij import WindowManager

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.gui.displaysettings.DisplaySettings import TrackMateObject
from fiji.plugin.trackmate.features.track import TrackIndexAnalyzer

import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer

from fiji.plugin.trackmate.io import TmXmlReader
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.io import CSVExporter
from fiji.plugin.trackmate.visualization.table import TrackTableView
from fiji.plugin.trackmate.action import ExportTracksToXML
from fiji.plugin.trackmate import Logger
from java.io import File

# ------------------------------
# 1) Read command-line arguments
#    argv[1] = full path to the stacked TIFF (e.g. /.../Stacked_tifs/Field0001_stack.tif)
#    argv[2] = path to output folder (e.g. /.../TrackMateResults/Field0001/)
# ------------------------------
if len(sys.argv) < 3:
    print("Usage: jython trackmate_jython_batch.py <stacked_tif_path> <output_folder>")
    sys.exit(1)

input_stack = sys.argv[1]
output_dir = sys.argv[2]

# Make sure output_dir exists
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

# Derive a “base name” from the TIFF:
# e.g. "/.../Stacked_tifs/Field0001_stack.tif" → base = "Field0001_stack"
base_name = os.path.splitext(os.path.basename(input_stack))[0]

# ------------------------------
# 2) Open the image in Fiji
# ------------------------------
imp = IJ.openImage(input_stack)
imp.show()

# ------------------------------
# 3) Create the TrackMate Model
# ------------------------------
model = Model()
model.setLogger(Logger.IJ_LOGGER)

# ------------------------------
# 4) Prepare Settings (same as your original) :contentReference[oaicite:4]{index=4}
# ------------------------------
settings = Settings(imp)

print("  nslices =", settings.nslices)
print("  nframes =", settings.nframes)
print("  zstart  =", settings.zstart)
print("  zend    =", settings.zend)
print("  tstart  =", settings.tstart)
print("  tend    =", settings.tend)

settings.nslices = 1
settings.nframes = 50
settings.zstart = 0
settings.zend = 0
settings.tstart = 0
settings.tend = 49

# Detector configuration:
settings.detectorFactory = LogDetectorFactory()
settings.detectorSettings = {
    'DO_SUBPIXEL_LOCALIZATION': True,
    'RADIUS': 0.025,
    'TARGET_CHANNEL': 1,
    'THRESHOLD': 4.0,
    'DO_MEDIAN_FILTERING': True,
}

# Spot filter (QUALITY > 0.7):
filter1 = FeatureFilter('QUALITY', 0.7, True)
settings.addSpotFilter(filter1)

# Tracker configuration:
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
settings.trackerSettings['ALLOW_TRACK_MERGING']  = False
settings.trackerSettings['ALLOW_GAP_CLOSING']   = True
settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 0.2
settings.trackerSettings['MAX_FRAME_GAP']       = 2
settings.trackerSettings['LINKING_MAX_DISTANCE']    = 0.2

# Add all analyzers:
settings.addAllAnalyzers()

# ------------------------------
# 5) Run TrackMate
# ------------------------------
trackmate = TrackMate(model, settings)
ok = trackmate.checkInput()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

ok = trackmate.process()
if not ok:
    sys.exit(str(trackmate.getErrorMessage()))

print("Process finished")

# ------------------------------
# 6) Display & log results (optional in headless)
# ------------------------------
selectionModel = SelectionModel(model)
ds = DisplaySettingsIO.readUserDefault()
ds.setTrackColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)
ds.setSpotColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)
displayer = HyperStackDisplayer(model, selectionModel, imp, ds)
#displayer.render()
displayer.refresh()
#model.getLogger().log(str(model))

# ------------------------------
# 7) Write XML + CSV into output_dir
# ------------------------------
# XML:
xml_filename = os.path.join(output_dir, base_name + ".xml")
xml_file = File(xml_filename)
writer = TmXmlWriter(xml_file, Logger.IJ_LOGGER)
writer.appendModel(model)
writer.appendSettings(settings)
writer.appendDisplaySettings(ds)
writer.appendGUIState("ConfigureViews")
writer.writeToFile()

# CSV (spots):
csv_spots = os.path.join(output_dir, base_name + ".csv")
CSVExporter.exportSpots(csv_spots, model, True)

# Spot table (only visible tracks):
spot_table = TrackTableView.createSpotTable(model, ds)
spot_table_csv = os.path.join(output_dir, base_name + "-spots.csv")
spot_table.exportToCsv(File(spot_table_csv))

# Edge table:
edge_table = TrackTableView.createEdgeTable(model, ds)
edge_table_csv = os.path.join(output_dir, base_name + "-edges.csv")
edge_table.exportToCsv(File(edge_table_csv))

# Track table:
track_table = TrackTableView.createTrackTable(model, ds)
track_table_csv = os.path.join(output_dir, base_name + "-tracks.csv")
track_table.exportToCsv(File(track_table_csv))

# Close the image:
imp.changes = False
imp.close()
