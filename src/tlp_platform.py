import pygmsh
import os
import tempfile
import numpy as np
from meshmagick import mesh as mm
from meshmagick import hydrostatics as hs
import capytaine as cpt
import logging
import xarray as xr
from scipy.linalg import block_diag


def pontoonArray(geom, draft, height, width, length, num):
    """
    Create a polar array of n pontoons

    Parameters
    ----------
    geom : TODO
        geometry engine from pygmsh
    draft : float/int
        vertical location of pontoons
    height : float/int
        pontoon height
    width : float/int
        pontoon width
    length : float/int
        pontoon length
    num : int
        number of pontoons

    Returns
    -------
    pontoons : TODO
        TODO

    """
    pontoons = list()
    for ii in range(0, num):
        pontoons.append(geom.add_box(x0=[-width/2, 0, draft+0*height/2],
                        extents=[width, length, height]))
        geom.rotate(pontoons[ii],
                    point=[0,0,0],
                    angle=ii/num*2*np.pi,
                    axis=[0,0,1])
    return pontoons


def calcImpedance(data, Bf=None):
    """
    Calculate hydro impedance

    Parameters
    ----------
    data : xarray.core.dataset.Dataset
        hydrodynamic and hydrostatic data from Capytaine
    Bf : TYPE, optional
        Frictional damping. The default is None.TODO

    Returns
    -------
    Zi : xarray.core.dataarray.DataArray
        impedance

    """
    if Bf is None:
        Bf = 0
    Zi = data.radiation_damping + \
            1j * (data.omega * (data.mass + data.added_mass) \
                  - data.hydrostatic_stiffness / data.omega )
    return Zi

class tlp_platform():
    """
    TLP platform object with methods for
    
    * generating geometry
    * running BEM to find hydrodynamic properties
    * getting hydrostatic properties    
    
    Parameters
    ----------
    r_spar : float/int
        spar radius
    draft : float/int
        platform draft
    height : float/int
        pontoon height
    width : float/int
        pontoon width
    length : float/int
        pontoon length
    num : int
        number of pontoons
    mass : float/int, optional
        mass of platform, if None, will be calculated based on displacement
    cog : list, optional
        location of center of gravity, if None, cog = [0, 0, -1*draft/2]
    mlambda : float/int, optional
        scaling factor, in None, mlambda = 1
    name : string, optional
        name for saved files, if None, name = 'platform_name'
    ofst : float/int, optional
        vertical offset, used to ensure that the mesh ends exactly at z=0, 
        if None, ofst = 0
    wrking_dir : string
        directory for saving files, if None, wrking_dir = '.' (current directory)
    """
    
    def __init__(self, r_spar, draft, height, width, length, num, mass=None, cog=None,
                 mlambda=None, name=None, ofst=None, wrking_dir=None):
        
        assert isinstance(r_spar, (int,float)), 'r_spar must be of type int or float'
        assert isinstance(draft, (int,float)), 'draft must be of type int or float'
        assert isinstance(height, (int,float)), 'height must be of type int or float'
        assert isinstance(width, (int,float)), 'width must be of type int or float'
        assert isinstance(num, int), 'num must be of type int'
        
        
        if mass is None:
            self.mass = None # TODO
        else:
            raise NotImplementedError() # TODO
            
        if cog is None:
            self.cog = [0, 0, -1*draft/2]
        else:
            raise NotImplementedError() # TODO    
            
        if mlambda is None:
            mlambda = 1
        
        if name is None:
            name = 'platform_name'

        if ofst is None:
            ofst = 0
            
        if wrking_dir is None:
            wrking_dir = '.'
        
        self.geom = pygmsh.opencascade.Geometry()
        
        self.name = name
        self.wrking_dir = wrking_dir
        
        self.ofst = ofst
        self.r_spar = r_spar
        self.draft = draft
        self.draft_ofst = self.ofst + draft
        
        self.height = height
        self.width = width
        self.length = length
        self.num = num
            
        spar = self.geom.add_cylinder(x0=[0.0, 0.0, self.ofst],
                                 axis=[0.0, 0.0, -1*self.draft_ofst],
                                 radius=self.r_spar,
                                 angle=2*np.pi,
                                 char_length=1)
        if num is not 0:
            pontoons = pontoonArray(geom=self.geom, 
                                    draft=-1*self.draft, 
                                    height=self.height, 
                                    width=self.width, 
                                    length=self.length, 
                                    num=self.num)
            
            # combine spar and pontoons
            self.geom.boolean_union(entities=pontoons + [spar], delete_first=True)
        
    def make_mesh(self, mshRefFactor=None, clcurv=None, write_file=None):
        if mshRefFactor is None:
            mshRefFactor = 1 # TODO
        
        if clcurv is None:
            clcurv = 360/50 # TODO 
            
        if write_file is None:
            write_file = True
            
            
        if write_file:
            out_dir = self.wrking_dir
        else:
            out_dir = dirpath = tempfile.mkdtemp()
            
        geo_filename = os.path.join(out_dir, self.name + '.geo')
        stl_filename = os.path.join(out_dir, self.name + '.stl')
        vtk_filename = os.path.join(out_dir, self.name + '.vtk')

            
        mshArgs = ['-clscale', str(mshRefFactor),
                   '-clcurv', str(clcurv)]
        
        self.mesh = pygmsh.generate_mesh(self.geom,
                                dim=2,
                                extra_gmsh_arguments=mshArgs,
                                remove_lower_dim_cells=True,
                                geo_filename=geo_filename,
                                msh_filename=stl_filename,
                                mesh_file_type="stl")
        
        self.mesh.write(vtk_filename)
        
    def run_hydro(self, dofs=None, freq=None, rho_water=None, wave_dirs=None):
        """
        Run calculations on platform to find hydrostatics and hydrodynamic
        properties

        Parameters
        ----------
        dofs : TODO, optional
            TODO The default is None.
        freq : np.array, optional
            frequencies for hydrodynamic calculations [Hz]. The default is np.linspace(5, 30, 5)
        rho_water : float, optional
            water density [kg/m3]. The default is 1025
        wave_dirs : TODO, optional
            wave directions for excitation calculations. The default is [0]

        Raises
        ------
        NotImplementedError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if freq is None:
            Tper = np.linspace(5, 30, 5)
            freq = 1/Tper
            
        if wave_dirs is None:
            wave_dirs = [0]

        omega = freq*2*np.pi
        
        if rho_water is None:
            self.rho_water = 1025
        else:
            self.rho_water = rho_water
        
        self.floatingbody = cpt.FloatingBody.from_file(self.name + '.stl', 
                                                       name=self.name)
        self.floatingbody.keep_immersed_part()
        
        if dofs is None:
            self.floatingbody.add_translation_dof(name='Surge')
            self.floatingbody.add_translation_dof(name='Heave')
            paxis = cpt.meshes.geometry.Axis((0, 1, 0), point=self.cog)
            self.floatingbody.add_rotation_dof(axis=paxis, name='Pitch')
        else:
            raise NotImplementedError() # TODO

            
        K, m, V, a = self.__getHydrostatics()
        self.floatingbody.mass = self.floatingbody.add_dofs_labels_to_matrix(block_diag(np.eye(2)*m,a.hs_data['Iyy']))
        self.floatingbody.hydrostatic_stiffness = self.floatingbody.add_dofs_labels_to_matrix(block_diag(0,np.array([[a.S33,a.S35],[a.S35, a.S55]])))
        
        
        self.hydro = self.__getHydrodynamics(omega=omega, wave_dirs=wave_dirs)
        self.hydro['displaced_volume'] = V
        
        
    def __getHydrostatics(self):
        """
        Find hydrostatic properties

        Returns
        -------
        K : TODO
            hydrostatic stiffness matrix
        m : TODO
            mass moment of inertia matrix
        V : float
            displaced volume
        a : meshmagick.Hydrostatics
            TODO

        """
        
        a = hs.Hydrostatics(mm.Mesh(self.floatingbody.mesh.vertices, self.floatingbody.mesh.faces),
                            cog=self.cog,
                            mass=self.mass,
                            rho_water=self.rho_water,
                            verbose=True)
        V = a.displacement_volume
        m = V * self.rho_water
        K = a.hydrostatic_stiffness_matrix
        return K, m, V, a
    
    
    def __getHydrodynamics(self, omega, wave_dirs, solver=None, log_level=None):
        """
        Run Capytaine boundary element method (BEM) solver to get hydrodynamic 
        properties.

        Parameters
        ----------
        omega : np.array
            wave frequencies [rad/s]
        wave_dirs : List[float]
            wave directions [TODO] 
        solver : TODO, optional
            BEM solver TODO The default is cpt.BEMSolver()
        log_level : TODO, optional
            TODO

        Raises
        ------
        NotImplementedError
            TODO

        Returns
        -------
        data : xarray.DataSet
            hydrodynamic data

        """
        
        if log_level is None:
            logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")
        
        test_matrix = xr.Dataset(coords={
            'omega': omega,
            'wave_direction': wave_dirs,
            'radiating_dof': list(self.floatingbody.dofs.keys()),
            })
        
        if solver is None:
            solver = cpt.BEMSolver()
        else:
            raise NotImplementedError() # TODO
            
        data = solver.fill_dataset(test_matrix, [self.floatingbody], 
                            hydrostatics=True, 
                            mesh=True, 
                            wavelength=True, 
                            wavenumber=True)
        
        data['freq'] = data.omega / (2 * np.pi)
        data['freq'].attrs['units'] = 'Hz'
        data = data.set_coords('freq')
        
        data['T'] = 1 / data.freq
        data['T'].attrs['units'] = 's'
        data = data.set_coords('T')
        
        data['Zi'] = calcImpedance(data)
        data['rao'] = cpt.post_pro.rao(data)

        return data
    
            
p = tlp_platform(r_spar=2, draft=30, height=5, width=2, length=10, num=0, ofst=1)
p.make_mesh(
    mshRefFactor=1,
    # clcurv=360/200, 
    write_file=True)
p.run_hydro()
