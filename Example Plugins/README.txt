EXAMPLE PLUGINS
===============

These files are REFERENCE ONLY. The app only loads plugins from:

    TMS EMG KIT/Fcns/Plugins/

To use a plugin, put a copy in the subfolder for its type. A plugin in the
wrong subfolder will not appear in the app.

    Fcns/Plugins/Custom_File_Import
        Reads a data file format the app doesn't support natively and converts
        it to the app's "Block" form (channels, blocks, sample rate, comments/
        events).

    Fcns/Plugins/Processing
        Include your own signal processing algorithms e.g. filtering,
        rectification, smoothing, detrending, artifact handling. 

    Fcns/Plugins/OnsetOffsetDetection
        Develop your own onset offset detection algorithms.
  

    Fcns/Plugins/Analysis
        Compute custom metrics from each trial within the determined onset and offset window
    
Each folder has a Template (start here) and simplified examples and various plugins in development based on manuscripts (will continue to update)

**Future development will include a way to share plugins with the community**


MAKING ONE
----------
1. Copy the Template for your plugin type and rename it. The function name
   inside must match the filename exactly - that's the name the app shows.
2. Edit only the marked ZONE blocks; leave the "DO NOT EDIT" section alone.
     ZONE 1  number of parameters, labels, defaults
     ZONE 2  unpack those values
     ZONE 3+ your method, filling the required outputs
3. Don't change the function signature. Return [] for outputs you don't compute.
4. Save it to the matching Fcns/Plugins subfolder and restart the app.

Each template's header lists the exact inputs, outputs, shapes and units.
