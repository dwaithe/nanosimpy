import numpy as np
import matplotlib.path as mplPath
import scipy.spatial as ssp
import scipy.stats as sst
import brownian_domain_trap_cy as bcy
import brownian_trap_cy as btcy
import copy

import equations_to_fit as fm
import sys







def calculate_psf(psfs, distance):
    """Calculates Gaussian of particular FWHM
    FWHM = 2*np.sqrt(2*np.log(2))*sigma                      # FWHM to sigma conversion
    sigma = FWHM/(2*np.sqrt(2*np.log(2)))                    # Sigma from FWHM
    G = np.exp(-x**2/(2*sigma**2))                           # is conventional Gaussian
    G = np.exp(-x**2/(2*(FWHM/(2*np.sqrt(2*np.log(2))))**2)) # substitute FWHM for sigma
    G = np.exp(-x**2/(2*(FWHM**2/(4*2*np.log(2)))))          # open the brackets and square contents
    G = np.exp((-x**2/(2*(FWHM**2)/8.))*(np.log(2)))         # decompose fraction
    G = np.exp((np.log(2.)))**(-x**2/((FWHM**2)/4.0))        # power law decomposition
    G = 2.**(-x**2/(FWHM/2.0)**2)                            # e^(ln2) = 2 indentity.
    
    """
 
    psf = {}
    psf['FWHMs'] = psfs
    psf['pixel_size'] = 1.0;
    psf['ri'] = np.meshgrid(np.arange(0,distance,psf['pixel_size']))[0];

    psf['number_FWHMs'] = psf['FWHMs'].__len__()
    psf['V'] = {}
    for ki in range(0, psf['number_FWHMs']):
        psf['V'][ki] = 2.0**(- psf['ri']**2 / (psf['FWHMs'][ki]/2.0)**2) #Gaussian function
    return psf
def integrate_over_psf_adv(psf,track_arr,num_of_mol,psy,psx,steps_to_jump,cc):
    #Pass each molecule through the psf function.

    
    psf['trace'] ={}
    sys.stdout.write('\n')
    for ki in range(0, psf['number_FWHMs']):
        sys.stdout.write('\r')
        sys.stdout.write("processing FWHM %d" % psf['FWHMs'][ki])
        sys.stdout.flush()
        trace = 0
        #sample track array at the desired time intervals.
        for b in range(0,num_of_mol):
            to_sample = np.round(np.arange(cc,track_arr[b][1].shape[0],steps_to_jump),0).astype(np.int32)
            track_x = track_arr[b][1][to_sample]
            track_y = track_arr[b][0][to_sample]
            trace += psf['V'][ki][np.round(np.sqrt((track_x-psx)**2+(track_y-psy)**2),0).astype(np.int32)]
        
        
        psf['trace'][ki] = copy.deepcopy(trace)
    return psf
def integrate_over_psf(psf,track_arr,num_of_mol,psy,psx):
    #Pass each molecule through the psf function.
    psf['trace'] ={}
    sys.stdout.write('\n')
    for ki in range(0, psf['number_FWHMs']):
        sys.stdout.write('\r')
        sys.stdout.write("processing FWHM %d" % psf['FWHMs'][ki])
        sys.stdout.flush()
        trace = 0
        for b in range(0,num_of_mol):
            trace += psf['V'][ki][np.round(np.sqrt((track_arr[b][1]-psx)**2+(track_arr[b][0]-psy)**2),0).astype(np.int32)]
        psf['trace'][ki] = copy.deepcopy(trace)
    return psf

def generate_hop_mesh(L, size, offset):
    seed_points =  int(np.round(4 * (offset / L)**2,0))
    points = (np.random.random((seed_points,2)))*size
    vor = ssp.Voronoi(points,qhull_options='Qbb Qc Qz')
    map_img = np.zeros((int(size),int(size)))
    area_arr = []
    vec_mast = []
    cc=1
    for reg in vor.regions:
        v_x = []
        v_y = []
        brek = False
        for b in reg:
            v0  = vor.vertices[b][0]
            v1  = vor.vertices[b][1]
            v_x.append(v0)
            v_y.append(v1)
            #print vor.vertices[b][1]
            if b == -1 or v0 <0 or v1 <0 or v0 >size or v1 >size:
                brek = True
        if brek == False and v_x !=[]:
            #Calculate area.
            
            area = 0
            for idx in range(2,v_x.__len__()):
                ptx_0 = v_x[0]
                pty_0 = v_y[0]
                ptx_1 = v_x[idx-1]
                pty_1 = v_y[idx-1]
                ptx_2 = v_x[idx]
                pty_2 = v_y[idx]
                
               
                a = np.sqrt((ptx_0-ptx_1)**2+(pty_0-pty_1)**2)
                b = np.sqrt((ptx_0-ptx_2)**2+(pty_0-pty_2)**2)
                c = np.sqrt((ptx_2-ptx_1)**2+(pty_2-pty_1)**2)
                p = (a + b + c)/2.0
                
                
                area += np.sqrt(p*(p-a)*(p-b)*(p-c))
            area_arr.append(area)
            
            vxmin = int(np.round(np.min(v_x),0))
            vxmax = int(np.round(np.max(v_x),0))
            vymin = int(np.round(np.min(v_y),0))
            vymax = int(np.round(np.max(v_y),0))
            ymesh, xmesh = np.meshgrid(np.arange(vymin,vymax),np.arange(vxmin,vxmax))
            dim = xmesh.shape
            bPath = mplPath.Path(np.array([v_y,v_x]).T)

            points_in =  bPath.contains_points(np.array([ymesh.reshape(-1),xmesh.reshape(-1)]).T)
            map_img[vymin:vymax,vxmin:vxmax] += points_in.reshape(dim).T*cc
            cc = cc+1
            v_x.append(v_x[0])
            v_y.append(v_y[0])
            vec_mast.append([v_x,v_y])
    print np.sqrt(np.average(area_arr))
    #Search the space for repetition.
    #narrow_list
    full_list = []
    for a in range(0,vec_mast.__len__()):
        vec = vec_mast[a]
        for b in range(1,vec[0].__len__()):
            full_list.append([vec[0][b-1],vec[0][b],vec[1][b-1],vec[1][b],0])
                    
            
    narrow_list = []
    dup_arr =[]
    for a in range(0,full_list.__len__()):
        vec = full_list[a]
        duplicate = False
        for b in range(0,narrow_list.__len__()):
            mec = narrow_list[b]
            
            if vec[0] == mec[0] and vec[1] == mec[1] and vec[2] == mec[2] and vec[3] == mec[3]:
                duplicate = True
                continue
                
            elif vec[0] == mec[1] and vec[1] == mec[0] and vec[2] == mec[3] and vec[3] == mec[2]:
                duplicate = True
                continue
            
                
        if duplicate == False:   
            narrow_list.append(vec)
    y1,x1 = np.where(map_img ==0)
    for x, y in zip(x1,y1):
        if x >100 and x <map_img.shape[1] - 100 and y>100 and y<map_img.shape[0]-100:
            
            med = np.max(map_img[y-5:y+5,x-5:x+5])
            
            if med > 0:
                map_img[y,x] = med
            
    return map_img, narrow_list, vec_mast, area_arr
def generate_trapping_zones(size, tzr,actual_number,num_positions_to_try,max_travel):
        
    #Divide this into the area.
    num_in_xy = np.ceil(np.sqrt(actual_number))

    #Create a meshgrid to distribute points equally across the full image size.
    ypoints, xpoints = np.meshgrid(np.arange(0,num_in_xy),np.arange(0,num_in_xy))
    ypoints = ypoints*(float(size)/float(num_in_xy))
    xpoints = xpoints*(float(size)/float(num_in_xy))
    ypoints = ypoints + float(tzr)
    xpoints = xpoints + float(tzr)
    #Round to nearest.
    ypoints = np.round(ypoints.reshape(-1),0).astype(np.int32)
    xpoints = np.round(xpoints.reshape(-1),0).astype(np.int32)

    #Convert to numpy array
    x_in = np.array(xpoints)
    y_in = np.array(ypoints)
    #We store locations as image map.
    loc_map = np.zeros((int(size),int(size)))
    #Possible direction combinatins
    pos_x = np.array([0, 0, 1, 1, 1,-1,-1,-1])
    pos_y = np.array([1,-1, 0, 1,-1, 0, 1,-1])
    #The empty positions are set as -1.
    loc_map = np.ones((size,size))*-1
    loc_map[y_in,x_in] = np.arange(0,x_in.shape[0]).astype(np.int32)

    

    #Our iteration counter.
    cc = 0
    
    num_of_iterations = 500

    one_per = int(np.round(num_of_iterations*actual_number*0.01, 0))
    total = float(num_of_iterations)*float(actual_number)

    while cc < total:
        
        if (cc % one_per) == 0.0:
                sys.stdout.write('\r')
                per = int(np.round(100*float(cc)/float(total),0))
                sys.stdout.write("Generating Map: [%-20s] %d%% complete" % ('='*(per/5), per))
                sys.stdout.flush()
        #Choose a random circle to move.
        cid = np.random.choice(int(actual_number),1)
        #Find the point.
        xc = x_in[cid]
        yc = y_in[cid]


        #Find and area around my coordinate.
        ycl = int(np.clip(yc-(tzr*4), 0, size))
        ycu = int(np.clip(yc+(tzr*4), 0, size))
        xcl = int(np.clip(xc-(tzr*4), 0, size))
        xcu = int(np.clip(xc+(tzr*4), 0, size))
        full_list = np.array(loc_map[ycl:ycu,xcl:xcu])
        
        #Generate a narrow list.
        narrow_list = full_list[full_list>-1]
        #If there are no neighbours within a range.
        if narrow_list.shape[0] == 0:
            #Number of iter
            for tti in range(0,num_positions_to_try):
                #Geneate random direction and distance.
                xco = pos_x[np.random.choice(np.arange(0,8),8,replace=False)]*np.random.choice(np.arange(0,int(max_travel)),1)
                yco = pos_y[np.random.choice(np.arange(0,8),8,replace=False)]*np.random.choice(np.arange(0,int(max_travel)),1)
                if x_in[cid] < size-1 and y_in[cid] < size-1 and x_in[cid]>1 and y_in[cid] >1:
                    
                    loc_map[y_in[cid], x_in[cid]] = -1
                    x_in[cid] += yco
                    y_in[cid] += xco
                    #Update the position on the location map.
                    loc_map[y_in[cid], x_in[cid]] = cid
                    break;
                else:
                    pass

        else:
            #The index cid is global lid
            b  = cid == narrow_list 
            x_in_cl = x_in[narrow_list.astype(np.int32)]
            y_in_cl = y_in[narrow_list.astype(np.int32)]
            #Copy input array.
            cli_x = list(copy.deepcopy(x_in_cl))
            cli_y = list(copy.deepcopy(y_in_cl))
            cli_xa = cli_x.pop(np.where(b)[0])
            cli_ya = cli_y.pop(np.where(b)[0])
            
            
            #Takes the first point that works. Skips if all iterations fail, circle didn't move.
            for tti in range(0,num_positions_to_try):
                #Random point and random direction.
                xco = pos_x[np.random.choice(np.arange(0,8),1)]*np.random.choice(np.arange(0,int(max_travel)),1)
                yco = pos_y[np.random.choice(np.arange(0,8),1)]*np.random.choice(np.arange(0,int(max_travel)),1)          
                
                #Try a combination
                cli_xb = cli_xa + xco
                cli_yb = cli_ya + yco
                
                #Make sure is a valid location.
                if cli_xb < size-1 and cli_yb < size-1 and cli_xb>1 and cli_yb >1:
                    
                    dx = cli_x- cli_xb
                    dy = cli_y- cli_yb
                    dist = np.sqrt((dx)**2+(dy)**2)
                    #Counts the number overlapping.
                    ctc = np.sum(dist < tzr*2)
                    
                    #If valid location.
                    if ctc == 0:
                        #Update the locations map.
                        loc_map[y_in[cid], x_in[cid]] = -1
                        loc_map[cli_yb,cli_xb] = cid
                        #Change the positions array
                        x_in[cid] = cli_xb
                        y_in[cid] = cli_yb
                        break

        cc += 1
    timg = np.zeros((int(size),int(size))).astype(np.int32)

    #Draw the spheres in their random locations.
    cc = 1

    total = x_in.shape[0]
    one_per = int(np.round(float(total)*0.01,0))
    sys.stdout.write('\n')
    for X, Y in zip(x_in.astype(np.int32),y_in.astype(np.int32)):
        if (cc % one_per) == 0.0:
                sys.stdout.write('\r')
                per = int(np.round(100*float(cc)/float(total),0))
                sys.stdout.write("Drawing Map: [%-20s] %d%% complete" % ('='*(per/5), per))
                sys.stdout.flush()
        #Search the radius around the object.
        for yc in range(-int(tzr*2),int(tzr*2)+1):
            for xc in range(-int(tzr*2),int(tzr*2)+1):
                    #If within one radius colour red.
                    if np.sqrt(yc**2 + xc**2) <= tzr:

                        timg[np.clip(yc+Y,0,int(size)-1),np.clip(xc+X,0,int(size)-1)] = cc
        cc += 1
    return timg, x_in, y_in


def brownian_spatial_hopping(total_sim_time,time_step,num_of_mol,D, map_img,R,offset,phob):
    D_out = D
    D_in = D
    phob_out = 1.0
    phob_in = phob
    return brownian_advanced_cython(total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out)
def brownian_stochastic_trap(total_sim_time,time_step,num_of_mol,D_out, D_in, R,offset,ptrap_off,ptrap_on):
    map_img = np.zeros((offset*2,offset*2))
    phob_out = 1.0
    phob_in = 1.0
    mode = 'brownian_trap'
    return brownian_advanced_cython(mode,total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out,ptrap_off,ptrap_on)
def brownian_stochastic_trap_hop(total_sim_time,time_step,num_of_mol,D_out, D_in,phob_in,phob_out,map_img, R,offset,ptrap_off,ptrap_on):
    
    mode = 'brownian_trap_hop'
    return brownian_advanced_cython(mode,total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out,ptrap_off,ptrap_on)
def brownian_only(total_sim_time,time_step,num_of_mol,D, R,offset):
    map_img = np.zeros((offset*2,offset*2))
    D_out = D
    D_in = D
    phob_out = 1.0
    phob_in = 1.0
    ptrap_off = -1
    ptrap_on =  -1
    mode = 'brownian'
    return brownian_advanced_cython(mode,total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out,ptrap_off,ptrap_on)

def brownian_only_numpy(total_sim_time,time_step,num_of_mol,D, width, height):
    """

    Inputs:
    total_simulation_time: Total simulation time in ms.
    time_step:             The duration of each time step ms.
    num_mol:               The number of molecules in the simulation.
    D:                     The diffusion rate.
    width:                 The width of the simulation.
    height:                The height of the simulation.
    Outputs:
    track_arr:             A dictionary where each track number (e.g. track_arr[0]) contains the track data [0,:] [1,:]
    """
   
    
    # Number of steps.
    num_of_steps = int(round(float(total_sim_time)/float(time_step),0))

    print 'num_of_steps',num_of_steps
    # Calculates length scales
    scale_in = np.sqrt(2.0 * (float(D)*1e3) * float(time_step))

    #Randomly generates start locations
    start_coord_x = (np.random.uniform(0.0, 1.0, num_of_mol))*width
    start_coord_y = (np.random.uniform(0.0, 1.0, num_of_mol))*height

    track_arr = {}
    #This can be done as one big matrix, but can crash system if large so 
    #I break it up by molecule.
    for b in range(0,num_of_mol):
        print 'processing tracks: ',(float(b)/float(num_of_mol))*100,'%'

        sys.stdout.write('\r')
        per = int((float(b)/float(num_of_mol))*100)
        sys.stdout.write("Processing tracks: [%-20s] %d%% complete" % ('='*(per/5), per))
        sys.stdout.flush()
        track = np.zeros((2,num_of_steps))
        track[0,0] = start_coord_y[b]
        track[1,0] = start_coord_x[b]
        rand_in  = sst.norm.rvs(size=[2,num_of_steps])*scale_in
        track[:,1:] += rand_in[:,1:]
        track = np.cumsum(track,1)
        out = track
        mod = np.zeros((out.shape))
        mod[0,:] = np.floor(track[0,:].astype(np.float64)/height)
        mod[1,:] = np.floor(track[1,:].astype(np.float64)/width)
        track_arr[b] = np.array(out-([mod[0,:]*height,mod[1,:]*width]))
        
        
    #We go through and make sure our particles wrap around.
    #for b in range(0,num_of_mol):
    #    print 'wrapping tracks: ',(float(b)/float(num_of_mol))*100,'%'
    #    bool_to_adapt = (track_arr[b][0,:]-offset)**2 + (track_arr[b][1,:]-offset)**2 >= R2
    #    while np.sum(bool_to_adapt) >0:
    #            ind = np.argmax(bool_to_adapt>0)
    #            phi = np.arctan2((track_arr[b][0,ind]-offset),(track_arr[b][1,ind]-offset));
    #            track_arr[b][1,ind:] = np.round(((track_arr[b][1,ind:]-offset) -  (2.0*(R-2) * np.cos(phi)))+offset,0).astype(np.int32);
    #            track_arr[b][0,ind:] = np.round(((track_arr[b][0,ind:]-offset) -  (2.0*(R-2) * np.sin(phi)))+offset,0).astype(np.int32);
    #            bool_to_adapt = (track_arr[b][0,:]-offset)**2 + (track_arr[b][1,:]-offset)**2 >= R2


    return track_arr
def brownian_spatial_trapping(total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out):
    mode = 'brownian_domain_trap'
    ptrap_off = -1
    ptrap_on =  -1
    return brownian_advanced_cython(mode,total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out,ptrap_off,ptrap_on)



def brownian_advanced_cython(mode,total_sim_time,time_step,num_of_mol,D_out,D_in, map_img,R,offset,phob_in,phob_out,ptrap_off,ptrap_on):
    """
    total_simulation_time: Total simulation time in ms.
    time_step:             The duration of each time step ms.
    num_mol:               The number of molecules in the simulation.
    D_in:                  The diffusion rate in diffusion zone.
    D_out:                 The diffusion rate outside the cell.
    Map_img:               An image with regions of hetrogeneity. background = 0 regions have index 1 to total_regions
    R:                     Radius of the simulation.
    offset:                Where the center of the simulation is occuring.
    phob_in:               Probability of hopping in to a region
    phob_out:              Probabliity of hopping out of a region.
    """

    # Number of steps.
    num_of_steps = int(round(float(total_sim_time)/float(time_step),0))
    # Calculates length scales
    scale_in = np.sqrt(2.0 * (float(D_in)*1e3) * float(time_step))
    scale_out = np.sqrt(2.0 * (float(D_out)*1e3) * float(time_step))
    

    #Randomly generates start locations
    t = np.random.uniform(0.0, 2.0*np.pi, num_of_mol)
    r = float(R) * np.sqrt(np.random.uniform(0.0, 1.0, num_of_mol))
    start_coord_x = r * np.cos(t) + float(offset)
    start_coord_y = r * np.sin(t) + float(offset)

    #Generate the region ids from the input map.
    map_img = map_img.astype(np.int32)
    id_list = map_img[(np.round(start_coord_y,0)).astype(np.int32),(np.round(start_coord_x,0).astype(np.int32))]
    trap_list = np.zeros((id_list.shape[0])).astype(np.int32)
    R2 = R**2
    width = map_img.shape[1]
    #We index as and array.
    map_img = (map_img.reshape(-1)).astype(np.int32)
    #This is the data structure that contains everything.
    track_arr ={}
    #Populate dictionary.
    for b in range(0,num_of_mol):
        track_arr[b] = np.zeros((2,num_of_steps))
        track_arr[b][0,0] = start_coord_y[b]
        track_arr[b][1,0] = start_coord_x[b]

    #For the loopy bit we use the cython.

    if mode == 'brownian_domain_trap':
        #Domain trapping needs its own dedicated function
        bcy.brownian_domain_trap(track_arr,num_of_steps, num_of_mol, id_list,width,map_img,phob_out,phob_in,R,R2,offset,scale_in,scale_out)
    
    if mode =='brownian_trap' or mode == 'brownian' or mode == 'brownian_trap_hop':
        #
    
        btcy.brownian_trap(track_arr,num_of_steps, num_of_mol,id_list,width,map_img,phob_out,phob_in,R,R2,offset,scale_in,scale_out,trap_list,ptrap_off,ptrap_on)
    return track_arr  