#!/usr/bin/python
# -*- coding: utf-8 -*-
##################################################
# © 2017 ETH Zurich, Swiss Seismological Service #
# Stefano Marano' - wavedec at gmail dot com   #
##################################################


import os, errno, sys
import yaml, csv
import logging
import argparse
import numpy as np
from numpy import pi, real, imag
from enum import Enum
import pulp as plp






class GRID(Enum):
  CARTESIAN = 1
  HEXAGONAL = 2
  RINGS_UNIFORM = 3
  RINGS_CUSTOM = 4
  RANDOM = 5
  FILE = 6



def placeSensors(Nsensors, pos, CostMask, normalized=False, EnforceSensor=None, TimeLimit=None, OUTPUT='.'):
  """
  The main function performing sensor placement
  
  Parameters:
    Nsensors        The number of sensors to be placed
    pos             Coordinate of the possible sensor positions
    CostMask        Where to sample the spatial frequency domain
    normalized      Normalize the cost function by the number of sensors
    EnforceSensors  list of position where to enforce the presence or absence of a snesor
    TimeLimit       time limit in seconds for the optimization
    OUTPUT          folder where to save any log file eg from solver
    
  Returns:
    pos_sol         Coordinate of the sensor positions
    status          Exit status of the solver

  """

  N = np.shape(pos)[0]
  M = np.shape(CostMask)[0]


  
  ###
  ### Define Fourier operator
  ###
  logging.debug("Creating Fourier operator")
  
  F = np.zeros((M,N)) + 1j*np.zeros((M,N))
  for mm in range(0,M):
    for nn in range(0,N):
      F[mm,nn] = np.exp(-2*pi*1j*np.dot(CostMask[mm,:], pos[nn,:]))
  if normalized:
    F /= Nsensors
  assert (M,N) == np.shape(F)
  
  pos_sol=None # solution to be returned if optimization fails
  try: 
    logging.debug("Setting up MIP model")
    SP_model = plp.LpProblem('SensorPlacement', plp.LpMinimize)
  
    # TODO    m.setParam("LogFile", os.path.join(OUTPUT,'./SP_gurobi.log'))
    
    # Create variables
    y = plp.LpVariable(cat=plp.LpContinuous, lowBound=0, name="y")
    x  = {(nn): plp.LpVariable(cat=plp.LpBinary, name="x_{0}".format(nn)) for nn in range(N)}
    SP_model.objective = plp.LpAffineExpression([ (y,1)])
  
  
    
    logging.debug("Setting up constraints")
    # Add inequality constraints
    for mm in range(M):
      SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum((+real(F[mm,nn]) * x[nn] for nn in range(N))- y), 
      sense=plp.LpConstraintLE, rhs=0, name="+real(F_{0})*x-y<=0".format(mm)))
  #      m.addConstr(LinExpr(+real(F[mm,:]), [x[nn] for nn in range(N)]), GRB.LESS_EQUAL, y, name="real(F)*x<=y")
    for mm in range(M):
      SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum((+imag(F[mm,nn]) * x[nn] for nn in range(N))- y), 
      sense=plp.LpConstraintLE, rhs=0, name="+imag(F_{0})*x-y<=0".format(mm)))
  #      m.addConstr(LinExpr(+imag(F[mm,:]), [x[nn] for nn in range(N)]), GRB.LESS_EQUAL, y, name="imag(F)*x<=y")
    for mm in range(M):
      SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum((real(F[mm,nn]) * x[nn] for nn in range(N))+ y), 
      sense=plp.LpConstraintGE, rhs=y, name="real(F_{0})*x+y>=0".format(mm)))
  ##      m.addConstr(LinExpr(-real(F[mm,:]), [x[nn] for nn in range(N)]), GRB.LESS_EQUAL, y, name="-real(F)*x<=y")
    for mm in range(M):
      SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum((imag(F[mm,nn]) * x[nn] for nn in range(N))+ y), 
      sense=plp.LpConstraintGE, rhs=y, name="imag(F_{0})*x+y>=0".format(mm)))
  ##      m.addConstr(LinExpr(-imag(F[mm,:]), [x[nn] for nn in range(N)]), GRB.LESS_EQUAL, y, name="-imag(F)*x<=y")
  #  
    SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum(x[nn] for nn in range(N)), 
      sense=plp.LpConstraintEQ, rhs=Nsensors, name="sum(x)=Nsensors".format(mm)))
  #    m.addConstr(LinExpr([1.0 for nn in range(N)], [x[nn] for nn in range(N)]), GRB.EQUAL, Nsensors, name="sum(x)=Nsensors")
    
    if EnforceSensor is not None:
      for nn in range(np.shape(EnforceSensor)[0]):
        d = np.linalg.norm(pos -  EnforceSensor[nn,0:2] , axis=1)
        ndx = np.argmin(d)
        d_min = d[ndx]
        if d_min < 1e3: # TODO set this intelligently
          SP_model.addConstraint(plp.LpConstraint(e=plp.lpSum(x[ndx]), 
            sense=plp.LpConstraintEQ, rhs=EnforceSensor[nn,2], name="EnforceSensor_{0}".format(nn)))
        else:
          logging.warning("Could not enforce sensor position ({},{})={}\nThe closest available position was {} away.".format(EnforceSensor[nn,0], EnforceSensor[nn,1], EnforceSensor[nn,2], d_min))
      
  
    EnforceMOI = False
    T = 5.1 # difficult to enforce exactly on a general grid
    # TODO need to enforce Qab < T or Qab = 0 for a couple of anlges.
    # split two cases for T == 0 and T > 0
    # need to introduce another variable to make the product of p_ij = x_i * x_j linear in the computaiton of Qab
    if EnforceMOI:
#      for psi in [0, np.pi/4]:
      expr_x = plp.lpSum(pos[nn,0] * x[nn] for nn in range(N))
      SP_model.addConstraint(plp.LpConstraint(e=expr_x, sense=plp.LpConstraintLE, rhs=+T, name="MOI_x_1"))
      SP_model.addConstraint(plp.LpConstraint(e=expr_x, sense=plp.LpConstraintGE, rhs=-T, name="MOI_x_2"))
      
      expr_y = plp.lpSum(pos[nn,1] * x[nn] for nn in range(N))
      SP_model.addConstraint(plp.LpConstraint(e=expr_y, sense=plp.LpConstraintLE, rhs=+T, name="MOI_y_1"))
      SP_model.addConstraint(plp.LpConstraint(e=expr_y, sense=plp.LpConstraintGE, rhs=-T, name="MOI_y_2"))
#          m.addConstr(LinExpr([(np.cos(psi)*pos[nn,0]+np.sin(psi)*pos[nn,1])*(-np.sin(psi)*pos[nn,0]+np.cos(psi)*pos[nn,1]) for nn in range(N)], [x[nn] for nn in range(N)]), GRB.EQUAL, 0.0, name="sum(x_n*y_n)=0")
#      print("MOI_x {}\n".format(pos[:,0]@x_sol))
#      print("MOI_y {}\n".format(pos[:,1]@x_sol))
#      print(pos_sol)
  
  
  
    logging.info("Initializing external solver")
    
    
    # solver = pulp.PULP_CBC_CMD(msg=1, mip_start=1)
    # solver = pulp.CPLEX_CMD(msg=1, mip_start=1)
    # solver = pulp.GUROBI_CMD(msg=1, mip_start=1)
    solver = plp.PULP_CBC_CMD(msg=True, mip=1, maxSeconds=TimeLimit) # TODO quite some options to set https://coin-or.github.io/pulp/technical/solvers.html#pulp.solvers.PULP_CBC_CMD
    
    logging.info("Running optimization")
    SP_model.solve(solver)
    status = plp.LpStatus[SP_model.status]
    logging.info("Optimization terminated. Status: {0}".format(status))
    
    if status == 'Optimal':
      #    y_sol = y.varValue
      x_sol = np.zeros((N,))
      for nn in range(N):
        x_sol[nn] = x[nn].varValue
  
      # solution
      ndx = np.where(x_sol > 0.5)[0]
      pos_sol = pos[ndx,:]
    
  except Exception as e:
    logging.critical('Error:' + str(e))
    status = plp.LpStatus[0]
    
  
  return pos_sol, status

def costMask(Kmin, Kmax, M):
  ###
  ### Define points sampled in frequency domain
  ###
  logging.debug("Defining sampling points in frequency domain")
  Area = pi*(Kmax**2 - Kmin**2)/2
  k_step = np.sqrt(Area / M)

  K_vec=np.linspace(Kmin, Kmax, np.int(np.max([np.floor((Kmax-Kmin)/k_step), 3])))

  CostMask_Azimuth = []
  CostMask_K = []
  for kk in K_vec:
    Azimuth_vec=np.linspace(0, pi, np.int(np.max([np.ceil(pi*kk/k_step), 4])), endpoint=False)
    for aa in Azimuth_vec:
      CostMask_Azimuth.append(aa)
      CostMask_K.append(kk)
  M=len(CostMask_Azimuth)
  CostMask_Azimuth = np.array(CostMask_Azimuth)
  CostMask_K = np.array(CostMask_K)
  
  CostMask = np.zeros((M,2))
  CostMask[:,0]=np.cos(CostMask_Azimuth)*CostMask_K
  CostMask[:,1]=np.sin(CostMask_Azimuth)*CostMask_K
  
  assert M == len(CostMask_Azimuth)
  assert M == len(CostMask_K)
  assert M == len(CostMask)
  
  return (CostMask, M)
  
def possiblePositions(N, L, Grid, GridFile=None, GridRings=None):

  ###
  ### Define arrangment of possible sensor positions
  ###
  logging.debug("Defining possible positions")
  if Grid == GRID.RANDOM: # random possible positions
    pos = np.random.rand(N,2)*L - L/2
  elif Grid == GRID.CARTESIAN: # Cartesian grid
    N = np.int(np.round(np.sqrt(N))**2)
    l_vec = np.linspace(-L/2, L/2, np.sqrt(N))
    pos_x, pos_y = np.meshgrid(l_vec, l_vec)
    pos = np.vstack((pos_x.flatten(), pos_y.flatten())).transpose()
  elif Grid == GRID.HEXAGONAL: # Hexagonal grid
    N = np.int(np.round(np.sqrt(N))**2)
    pos = np.zeros((N,2))
    l_step = L/(np.sqrt(N)-1)
    l_vec = np.linspace(-L/2, L/2, np.sqrt(N))
    nn = 0
    for xx in range(np.int(np.sqrt(N))):
      for yy in range(np.int(np.sqrt(N))):
        pos[nn,0] = -L/2 + xx*l_step
        pos[nn,1] = -L/2 + yy*l_step + np.mod(xx,2)*l_step/2
        nn +=1
  elif Grid == GRID.RINGS_UNIFORM:
    pos = []
    Radius = L/2
    Area = pi*Radius**2
    r_step = np.sqrt(Area / N)
  
    R_vec=np.linspace(0, Radius, np.int(np.max([np.ceil(Radius/r_step), 2])))

    for rr in R_vec:
      Azimuth_vec=np.linspace(0, 2*pi, np.int(np.max([np.ceil(2*pi*rr/r_step), 4])), endpoint=False)
      if rr == 0: Azimuth_vec= np.array([0])
      for aa in Azimuth_vec:
        pos.append(np.array([rr*np.cos(aa), rr*np.sin(aa)]))
    pos = np.array(pos)
    N = len(pos)
  elif Grid == GRID.RINGS_CUSTOM: # Circular rings
    if GridRings is None:
      logging.critical("Parameter GridRings not found")
      return
    pos = []
    Radius = L/2

    Grid_r = Radius*GridRings[:,0]  # GridRings[:,0] contains r1,r2,r3... radii normalized between 0 and 1
    Grid_n = GridRings[:,1]     # GridRings[:,1] contains n1,n2,n3... number of sensors
    if GridRings[0,0] == 0:
      Delta_r = np.diff(Grid_r)
      r_step = np.sum(Delta_r * Grid_n[1]) / N
    else:
      Delta_r = np.concatenate(([0],Grid_r))
      r_step = np.sum(Delta_r * Grid_n) / N
      
    
    R_vec=np.linspace(0, Radius, np.int(np.max([np.ceil(Radius/r_step), 2])))
    
    # TODO anything to debug here?
#    print(Grid_r)
#    print(Grid_n)
    
    for radius in R_vec:
      
      n = next(n for n,r in enumerate(Grid_r) if r >= radius)
      Nr = Grid_n[n]
#      print("")
#      print(radius)
#      print(Nr)s
      if radius == 0:
        Nr = np.min([Nr, 1])
      Nr = np.int(Nr)
      for nr in range(0,Nr):
        pos.append(np.array([radius*np.cos(2*pi*nr/Nr), radius*np.sin(2*pi*nr/Nr)]))
    pos = np.array(pos)
    N = len(pos)
    
  elif Grid == GRID.FILE:
    if not os.path.exists(GridFile):
      logging.critical("File not found - GridFile: {0}.".format(GridFile))
      return
    with open(GridFile, 'r', encoding='UTF8') as f:
      reader = csv.reader(f, delimiter='\t')
      pos = []
      for row in reader:
        if not row[0].startswith('#') and len(row) == 2:
          pos.append(np.array([np.float(row[0]), np.float(row[1])]))
    pos = np.array(pos)
    N = len(pos)
  else:
    pos = None
    logging.critical("Grid type not chosen")
    return
  
  assert (N,2) == np.shape(pos)

  return (pos, N)

def main():

  

  ### Initialization: commandline args, config file, and other parameters
  parser = argparse.ArgumentParser(description='SP - design of seismic sensor array')
  parser.add_argument("--config_file", default='config.yaml', help="path to YAML configuration file [config.yaml]")
  parser.add_argument("--output", help="path to output folder [./]")
#  parser.add_argument("--Nsensors", help="Parameter Nsensors")
#  parser.add_argument("--Kmin", help="Parameter Kmin")
#  parser.add_argument("--Kmax", help="Parameter Kmax")
  parser.add_argument("--verbosity", type = int, default=1, help="increase output verbosity: 0 = only warnings, 1 = info, 2 = debug.  [0]")
  args = parser.parse_args()

  # Setup console logging, file logging later
  logging.basicConfig(level=logging.INFO,  format='%(message)s', stream=sys.stderr)

  
  # parameter CONIFG_FILE (configuration file in YAML format)
  CONFIG_FILE = os.path.abspath(args.config_file)  # default is "./config.yaml"
  try:
    with open(CONFIG_FILE) as f:
      conf = yaml.load(f, Loader=yaml.FullLoader)
      f.close()
  except yaml.scanner.ScannerError as e:
    logging.critical('There is a syntax problem with the YAML configuration file {0}.'.format(CONFIG_FILE))
    logging.critical(e)
    return
  except FileNotFoundError:
    logging.warning('No configuration file ({0}) found. Proceeding with default values.'.format(CONFIG_FILE))
    conf = dict()
  except Exception as e:
    logging.critical('There was a problem while loading the YAML configuration file {0}.'.format(CONFIG_FILE))
    logging.critical(type(e))
    logging.critical(e)
    return

  if conf is None: # when we read an empty config file
    conf = dict()
    
    
  # parameter OUTPUT  
  OUTPUT_args = args.output # the value of OUTPUT as parsed from commandline
  OUTPUT_conf = conf.get('OUTPUT', None) # the value of OUTPUT as read from CONFIG_FILE
  if (OUTPUT_args == None) and (OUTPUT_conf == None):    OUTPUT = '.'
  elif (OUTPUT_args == None) and (OUTPUT_conf != None):  OUTPUT = OUTPUT_conf
  elif (OUTPUT_args != None) and (OUTPUT_conf == None):  OUTPUT = OUTPUT_args
  else:
    logging.critical('Option clash for \'OUTPUT\'. Defined both as commandline argument and in config file \'{0}\''.format(CONFIG_FILE))
    return
  OUTPUT = os.path.abspath(OUTPUT)
  
  # create output folder and make sure it is writable
  try:
    os.makedirs(OUTPUT)
  except OSError as exception:
    if exception.errno != errno.EEXIST:  raise
  
  
  # logging to file    
  logFile = logging.FileHandler(filename=os.path.join(OUTPUT, 'SP.log'), mode='w')
  logFile.setLevel(logging.DEBUG)
  formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', "%Y-%m-%d %H:%M:%S")
  logFile.setFormatter(formatter) # tell the handler to use this format
  logging.getLogger('').addHandler(logFile) # add the handler to the root logger
  
  logging.debug('Starting SP with Python version:')
  logging.debug(sys.version_info)
  logging.debug('NumPy version: {0}'.format(np.__version__))
  
  # parsing YAML configuration file  
  try:
    # Mandatory parameters
    Nsensors = conf.get('Nsensors', None)
    Kmin = conf.get('Kmin', None)
    Kmax = conf.get('Kmax', None)
    if Nsensors is None:
      logging.critical("Parameter Nsensors not found in configuration file {0}".format(CONFIG_FILE))
      return
    if Kmin is None:
      logging.critical("Parameter Kmin not found in configuration file {0}".format(CONFIG_FILE))
      return
    if Kmax is None:
      logging.critical("Parameter Kmax not found in configuration file {0}".format(CONFIG_FILE))
      return
    Nsensors = int(Nsensors)
    Kmin = float(Kmin)
    Kmax = float(Kmax)
    if Nsensors < 3:
      logging.critical("Parameter Nsensors must be >= 3.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
    if Kmin <= 0.0:
      logging.critical("Parameter Kmin must be > 0.0.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
    if Kmax <= Kmin:
      logging.critical("Parameter error, must be Kmax >= Kmin.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
      
    # Optional parameters
    L = conf.get('MaximumAperture', None)
    if L is None:
      L = 1.2/Kmin
    else:
      L = np.float(L)
    if L <= 0.0:
      logging.critical("Parameter L must be > 0.0.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
    N = np.int(conf.get('N', 200)) # this is intended, not actual
    M = np.int(conf.get('M', 100))  # this is intended, not actual
    if N <= Nsensors:
      logging.critical("Parameter N too small.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
    if M <= 0:
      logging.critical("Parameter M must be > 0.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return

    GridFile = conf.get('GridFile', None)
    GridRings = conf.get('GridRings', None)
    Grid = conf.get('Grid', None)
    
    if GridFile is not None and str(Grid).lower() != 'file':
      logging.critical("Parameter GridFile specified but Grid is not set to 'file'.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
    if GridRings is not None and str(Grid).lower() != 'rings_custom':
      logging.critical("Parameter GridRings specified but Grid is not set to 'rings_custom'.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
     
    Grid = str(Grid).lower()
    if Grid == 'cartesian':
      Grid = GRID.CARTESIAN
    elif Grid == 'hexagonal':
      Grid = GRID.HEXAGONAL
    elif Grid == 'rings_uniform':
      Grid = GRID.RINGS_UNIFORM
    elif Grid == 'rings_custom':
      Grid = GRID.RINGS_CUSTOM
    elif Grid == 'random':
      Grid = GRID.RANDOM
    elif Grid == 'file':
      Grid = GRID.FILE
    elif Grid == str(None).lower():
      Grid = GRID.RINGS_UNIFORM
    else:
      logging.critical("Error while reading parameter Grid.")
      logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
      return
      
    
    if Grid == GRID.FILE:
      GridFile = str(GridFile)
      if not os.path.exists(GridFile):
        logging.critical("File not found - GridFile: {0}.".format(GridFile))
        logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
        return
    elif Grid == GRID.RINGS_CUSTOM:
      if GridRings is None:
        logging.critical("Please supply parameter GridRings. It must be supplied when Grid: 'rings_custom'")
        logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
        return
      GridRings = np.array(GridRings)
      if len(np.shape(GridRings)) != 2:
        logging.critical("Check parameter GridRings. Must be a two dimensional array [[r1,n1],[r2,n2],[r3,n3],...]")
        logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
        return
      # TODO check r is between 0 and 1 (first element must be 0, last element must be 1), check n is integer

    
    EnforceSensor = conf.get('EnforceSensor', None) # 0/1 constrained else unconstrained
    if EnforceSensor is not None:
      EnforceSensor = np.array(EnforceSensor)
      if len(np.shape(EnforceSensor)) != 2:
        logging.critical("Check parameter EnforceSensor. Must be a two dimensional array [[x1,y1,0],[x2,y2,1],[x3,y3,0],...]")
        logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
        return
      if not np.all(np.logical_or(EnforceSensor[:,2]  == 0.0 , EnforceSensor[:,2]  == 1.0)):   
        logging.critical("Check parameter EnforceSensor. Third column must be 0/1")
        logging.critical("Check configuration file {0}.".format(CONFIG_FILE))
        return
      
    TimeLimit = np.int(conf.get('TimeLimit', 120.0))
    
    normalized = False
    

    
  except Exception as e:
    logging.critical('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))
    logging.critical(e)
    raise ValueError('There was a problem while parsing the configuration file {0}'.format(CONFIG_FILE))




  logging.info("")
  logging.info("Sensor Placement for seismic surface waves")
  logging.info("© 2017 ETH Zurich, Swiss Seismological Service")
  logging.info("")
  logging.info("\tWith configuration file")
  logging.info("\t\t{0}".format(CONFIG_FILE))
  logging.info("\tOutput will be stored in")
  logging.info("\t\t{0}".format(OUTPUT))
 

  (pos, N) = possiblePositions(N, L, Grid, GridFile, GridRings)
  (CostMask, M) = costMask(Kmin, Kmax, M)
  
  logging.info("\tPlacing {0} sensors".format(Nsensors))
  logging.info("\tWith {0} possible sensor positions".format(N))
  logging.info("\tMaximum aperture {0:.2e}  [m]".format(L))
  logging.info("\tKmin={0:.2e} Kmax={1:.2e} [1/m]".format(Kmin, Kmax))
  
  with open(os.path.join(OUTPUT, 'PossiblePositions.csv'), 'w', encoding='UTF8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['# Coordinates of possible sensor positions'])
    writer.writerow(['# N={0}'.format(N)])
    writer.writerow(['# Easting', 'Northing'])
    writer.writerow(['# [m]', '[m]'])
    for nn in range(N):
      writer.writerow(['{:6.3e}'.format(pos[nn,0]), '{:6.3e}'.format(pos[nn,1])])

  with open(os.path.join(OUTPUT, 'CostMask.csv'), 'w', encoding='UTF8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['# Locations of frequency sampling points'])
    writer.writerow(['# M={0}'.format(M)])
    writer.writerow(['# Wavenumber x', 'Wavenumber y'])
    writer.writerow(['# [1/m]', '[1/m]'])
    for mm in range(M):
      writer.writerow(['{:6.3e}'.format(CostMask[mm,0]), '{:6.3e}'.format(CostMask[mm,1])])


  
  # TODO maybe a loop would be nice, but there are other priorities
  pos_sol, status = placeSensors(Nsensors, pos, CostMask, normalized=normalized, EnforceSensor=EnforceSensor, TimeLimit=TimeLimit, OUTPUT=OUTPUT)
  

  if pos_sol is None:
    logging.critical("Optimization failed. No array was generated.")
    return
  else:
    with open(os.path.join(OUTPUT, 'Optimized_ArrayLayout.csv'), 'w', encoding='UTF8') as f:
      writer = csv.writer(f, delimiter='\t')
      writer.writerow(['# Coordinates of optimized array'])
      writer.writerow(['# Easting', 'Northing'])
      writer.writerow(['# [m]', '[m]'])
      for nn in range(Nsensors):
        writer.writerow(['{:6.3e}'.format(pos_sol[nn,0]), '{:6.3e}'.format(pos_sol[nn,1])])
    d = {'Nsensors' : Nsensors,
       'Kmin' : Kmin, 'Kmax' : Kmax,
       'N' : N, 'M': M,
       'MaximumAperture' : L, 'TimeLimit' : TimeLimit}
    # if EnforceSensor is not None:
    #  d['EnforceSensor'] = EnforceSensor
    if Grid == GRID.RANDOM: # random possible positions
      d['Grid'] = 'random'
    elif Grid == GRID.CARTESIAN: # Cartesian grid
      d['Grid'] = 'cartesian'
    elif Grid == GRID.HEXAGONAL: # Hexagonal grid
      d['Grid'] = 'hexagonal'
    elif Grid == GRID.RINGS_UNIFORM:
      d['Grid'] = 'rings_uniform'
    elif Grid == GRID.RINGS_CUSTOM: # Circular rings
      d['Grid'] = 'rings_custom'
    elif Grid == GRID.FILE:
      d['Grid'] = 'file'
      d['GridFile'] = GridFile
      
    with open(os.path.join(OUTPUT, 'Optimized_Info.yaml'), 'w') as outfile:
      yaml.dump(d, outfile, default_flow_style=False)
  
    if True: # TODO remove or keep this preliminary plot?
      import matplotlib.pyplot as plt
      fig = plt.figure()
      plt.plot(pos[:,0], pos[:,1], 'o', ms=2, lw=0, alpha=1.0, mfc='blue')
      plt.plot(pos_sol[:,0], pos_sol[:,1], 'o', ms=4, lw=0, alpha=1.0, mfc='red')
      ax = fig.gca()
      ax.set_aspect('equal')
      ax.set_xlabel('Easting [m]')
      ax.set_ylabel('Northing [m]')
      xMin = np.min(pos[:,0]); xMax = np.max(pos[:,0]);
      yMin = np.min(pos[:,1]); yMax = np.max(pos[:,1]);
      xAperture = np.max([xMax - xMin, 0.8*(yMax - yMin)]) # avoid really stretched aspect ratios
      yAperture = np.max([yMax - yMin, 0.8*(xMax - xMin)])
      ax.set_xlim(xMin-0.05*xAperture, xMax+0.05*xAperture)
      ax.set_ylim(yMin-0.05*yAperture, yMax+0.05*yAperture)
      plt.pause(15)
  
  return



if  __name__ =='__main__':
  main()

