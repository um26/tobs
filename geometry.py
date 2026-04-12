from ngsolve import Mesh
from netgen.geom2d import SplineGeometry

def transformer(n_coil : int = 2,
                # dimensions
                width_box : float = 0.6,
                height_box : float = 0.6,
                inner_diameter_coil : float = 4e-2,
                outer_diameter_coil : float = 6e-2,
                height_coil : float = 4e-2,
                # labels
                primary_plus : str = "Pp",
                primary_minus : str = "Pm",
                secondary_plus : str = "Sp",
                secondary_minus : str = "Sm",
                design_domain : str = "Omega_c",
                bc_out : str = "dOmega",
                bc_coil_primary_plus : str = "dPp",
                bc_coil_primary_minus : str = "dPm",
                bc_coil_secondary_plus : str = "dSp",
                bc_coil_secondary_minus : str = "dSm",
                # mesh
                maxh : float = 1e-2,
                maxh_coil : float = 1e-2,
                maxh_design : float = 1e-2
                ) -> Mesh:
    """ Generates the mesh of the transformer geometry """

    ## Geometry generation
    geo = SplineGeometry()
    geo.AddRectangle( p1 = (-width_box/2, -height_box/2),
                      p2 = ( width_box/2,  height_box/2),
                      leftdomain = 1,
                      rightdomain = 0,
                      bc = bc_out)

    if n_coil == 1:
        geo.AddRectangle( p1 = (-outer_diameter_coil/2, -height_coil/2), p2 = (-inner_diameter_coil/2, height_coil/2),
        leftdomain = 2, rightdomain = 1, bc = bc_coil_primary_plus )
        geo.AddRectangle( p1 = (inner_diameter_coil/2, -height_coil/2), p2 = (outer_diameter_coil/2, height_coil/2),
        leftdomain = 3, rightdomain = 1, bc = bc_coil_primary_minus)
        # Set materials labels
        geo.SetMaterial( 2, primary_plus)      # Primary conductor (positive)
        geo.SetMaterial( 3, primary_minus)     # Primary conductor (negative)
        # Set mesh sizes
        geo.SetDomainMaxH( 2, maxh_coil )    
        geo.SetDomainMaxH( 3, maxh_coil ) 

    elif n_coil == 2:
        th = (outer_diameter_coil - inner_diameter_coil) / 2
        ep = (outer_diameter_coil + inner_diameter_coil) / 2
        geo.AddRectangle( p1 = (-width_box/2 + (0.3-0.11), -height_coil/2), p2 = (-width_box/2 + (0.3-0.11) + th, height_coil/2),
        leftdomain = 2, rightdomain = 1, bc=bc_coil_primary_plus )
        geo.AddRectangle( p1 = (-width_box/2 + (0.3-0.11) + ep, -height_coil/2), p2 = (-width_box/2 + (0.3-0.11) + ep + th, height_coil/2),
        leftdomain = 3, rightdomain = 1, bc=bc_coil_primary_minus )
        geo.AddRectangle( p1 = (width_box/2 - (0.3-0.11) - th - ep, -height_coil/2), p2 = (width_box/2 - (0.3-0.11) - ep, height_coil/2),
        leftdomain = 4, rightdomain = 1, bc=bc_coil_secondary_minus )
        geo.AddRectangle( p1 = (width_box/2 - (0.3-0.11) - th, -height_coil/2), p2 = (width_box/2 - (0.3-0.11) , height_coil/2),
        leftdomain = 5, rightdomain = 1, bc=bc_coil_secondary_plus )
        # Set materials labels
        geo.SetMaterial( 2, primary_plus )      # Primary conductor (positive)
        geo.SetMaterial( 3, primary_minus )     # Primary conductor (negative)
        geo.SetMaterial( 4, secondary_minus )   # Secondary conductor (negative)
        geo.SetMaterial( 5, secondary_plus )    # Secondary conductor (positive)
        # Set mesh sizes
        geo.SetDomainMaxH( 2, maxh_coil )    
        geo.SetDomainMaxH( 3, maxh_coil )  
        geo.SetDomainMaxH( 4, maxh_coil )   
        geo.SetDomainMaxH( 5, maxh_coil )

    # Name the design domain and set its mesh size
    geo.SetMaterial( 1, design_domain ) 
    geo.SetDomainMaxH(1, maxh_design )

    return Mesh( geo.GenerateMesh( maxh = maxh ))