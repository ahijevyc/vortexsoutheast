import glob
from metpy.units import units
import os
import spc
import sys

workdir = "/glade/scratch/ahijevyc/vortexsoutheast/stormtimelist"

def composites():
    # get all files ending with .00z.txt
    ff = glob.glob(os.path.join(workdir,"*z.txt"))
    # split filename by periods and get first part
    cl = [os.path.basename(f).split('.')[0] for f in ff]
    return set(cl)

def str3(label,az,r):
    # stringify 3 values with slash delimiter
    return f"{label}/{az}deg/{r}km"

def ptlist(label,az,r):
    # return point label and coordinates at centroid and opposite centroid
    opposite_az = (az + 180) % 360
    return [str3(label,az,r), str3("opposite "+label, opposite_az, r)]

# dictionary with coords as keys and a dictionary of pointlists as values.
pointlist = {
        # Ran NARR_composite.py with all LTC categories except near coast and well inland,
        # 0003z, 0609z, 1215z, 1821z combined (/glade/work/ahijevyc/share/VSE/all_LTC.txt)
        # these are points in the composite corresponding to max shear or torn rpts
        "north" : {
                "shr10_700 max"                  : ptlist("shear max", 61.5, 260), 
                "torn max"                       : ptlist("torn max",  74.3, 254), 
                "tornadoes_well_inland torn max" : ptlist("torn max",  80,   200)
                },
        "storm motion" : {
                "shr10_700 max"                  : ptlist("shear max", 47.5, 260), 
                "torn max"                       : ptlist("torn max",  60.0, 245), 
                },
        "wind shear" : {
                "shr10_700 max"                  : ptlist("shear max",  2.5, 260), 
                "torn max"                       : ptlist("torn max",  19.2, 239), 
                },

        }


def centroids():
    # Return all possible centroid strings in pointlist dictionary.
    c = []
    for coord in pointlist:
        c.extend(pointlist[coord].keys()) 
    return set(c)

def tc_coords(df, trackdf):
    """
    get distance from TC center for each storm report in df
    """
    i = trackdf.valid_time == df.name
    assert i.sum(), f"no track times match {df.name}"
    originlon = trackdf.loc[i,"lon"].values[0] * units.degrees_E
    originlat = trackdf.loc[i,"lat"].values[0] * units.degrees_N
    dist_from_origin, heading = spc.gdist_bearing(originlon, originlat, 
            df["slon"].values * units.degrees_E, df["slat"].values * units.degrees_N)
    df["dist_from_origin"] = dist_from_origin
    return df


