def get_fiducial_info(raw_info, save_name):
    """Save fiducials and accompanying info as mat file.


    Parameters
    ----------
    raw_info : MNE-Python info object
        info of raw object, containing all relevant information
    save_name : str
        name or path to save the resulting MATLAB structure to

    Returns
    -------
    Saves a MATLAB file containing the following MATLAB objects:
    dev_head_t : MATLAB structure
        contains the device to head transform matrix
    fid_coords : MATLAB matrix (double)
        contains the coordinates of the fiducials in device space
    fid_id : MATLAB cell
        contains the fiducial identifiers for the rows of fid_coords

    Note: Tested only for Elekta data so far.
    """
    import scipy.io as spo

    dig_info = raw_info['hpi_results'][0]['dig_points']

    ident = []
    coords = np.zeros((len(dig_info), 3))
    for ii, points in enumerate(dig_info):
        ident.append(points['ident'])
        coords[ii, :] = points['r']

    mapping = {1 : 'lpa', 2 : 'nas', 3 : 'rpa', 4: 'unknown'}
    fid_id = np.zeros((len(ident), ), dtype=np.object)
    fid_id[:] = [mapping[id] if id in mapping else 'unknown' for id in ident]


    fiducial_info = {'fid_id': fid_id, 'fid_coord': coords,
                     'dev_head_t': raw.info['hpi_results'][0]['coord_trans']}

    spo.savemat(save_name, fiducial_info)
