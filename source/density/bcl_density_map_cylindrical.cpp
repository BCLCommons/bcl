// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "density/bcl_density_map_cylindrical.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_cylinder_coordinates.h"
#include "density/bcl_density_map.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MapCylindrical::s_Instance
    (
      GetObjectInstances().AddInstance( new MapCylindrical())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    MapCylindrical::MapCylindrical() :
      m_Body(),
      m_HeightResolution(),
      m_RadiusResolution(),
      m_AngleResolution(),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
    }

    //! construct MapCylindrical from given Parameters
    //! @param BODY main axis of cylindrical density map as linesegment3D
    //! @param DENSITY_MAP the original euclidean density map from which the cylindrical density map is constructed
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    MapCylindrical::MapCylindrical
    (
      const assemble::SSEGeometryInterface &BODY,
      const Map &DENSITY_MAP,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    ) :
      m_Body( BODY),
      m_HeightResolution( HEIGHT_RESOLUTION),
      m_RadiusResolution( RADIUS_RESOLUTION),
      m_AngleResolution( 2 * math::g_Pi / NUMBER_WEDGES),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
      BCL_Assert( BODY.GetOrientation().IsDefined(), "Body has to be defined");
      BCL_Assert
      (
        HEIGHT_RESOLUTION > double( 0) && RADIUS_RESOLUTION > double( 0) && NUMBER_WEDGES > 0,
        "height and radius resolution have to be larger than 0, number of wedges has to be larger than 0"
      );
      BCL_Assert( UPPER_RADIUS > RADIUS_RESOLUTION, "maximal radius has to be greater than radius resolution");

      //convert parameter density map in TriCubicSpline
      math::TricubicSpline densitymap_as_spline( DENSITY_MAP.ConvertDensityToSpline());

      //calculate the number of voxel in height, radius and angle
      const size_t number_height_voxels( size_t( 2 * BODY.GetExtent( coord::GetAxes().e_Z) / HEIGHT_RESOLUTION));
      const size_t number_radius_voxels( size_t( UPPER_RADIUS / RADIUS_RESOLUTION));
      const size_t number_angle_voxels( NUMBER_WEDGES);

      // resize data (create tensor of appropriate size, filled with zeros)
      m_Data =
        math::Tensor< double>
       (
         number_height_voxels,
         number_radius_voxels,
         number_angle_voxels,
         double( 0)
       );

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < number_height_voxels; ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( 0); itr_radius < number_radius_voxels; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < number_angle_voxels; ++itr_angle)
          {
            // create cylindrical coordinates
            const coord::CylinderCoordinates cylind_coord
            (
              ( int( itr_length) - int( number_height_voxels / 2)) * m_HeightResolution,
              itr_radius * m_RadiusResolution,
              itr_angle * m_AngleResolution
            );

            // convert to Cartesian coordinates
            linal::Vector3D cart_coord( cylind_coord.GetCartesianCoordinates());

            // transform Cartesian coordinates with Body
            cart_coord.Transform( BODY.GetOrientation());

            // lookup intensity in spline and insert into this density map
            m_Data( itr_length, itr_radius, itr_angle) =
              densitymap_as_spline.F( cart_coord.Z(), cart_coord.Y(), cart_coord.X());
          }
        }
      }
    }

    //! construct MapCylindrical from given Parameters
    //! @param BODY
    //! @param SPLINE the spline of the original euclidean density map
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    MapCylindrical::MapCylindrical
    (
      const assemble::SSEGeometryInterface &BODY,
      math::TricubicSpline &SPLINE,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    ) :
      m_Body( BODY),
      m_HeightResolution( HEIGHT_RESOLUTION),
      m_RadiusResolution( RADIUS_RESOLUTION),
      m_AngleResolution( 2 * math::g_Pi / NUMBER_WEDGES),
      m_Data(),
      m_Minimum(),
      m_Maximum(),
      m_Mean()
    {
      BCL_Assert( BODY.GetOrientation().IsDefined(), "Body has to be defined");
      BCL_Assert
      (
        HEIGHT_RESOLUTION > double( 0) && RADIUS_RESOLUTION > double( 0) && NUMBER_WEDGES > 0,
        "height and radius resolution have to be larger than 0, number of wedges has to be larger than 0"
      );
      BCL_Assert( UPPER_RADIUS > RADIUS_RESOLUTION, "maximal radius has to be greater than radius resolution");

      //calculate the number of voxel in height, radius and angle
      const size_t number_height_voxels( size_t( 2 * BODY.GetExtent( coord::GetAxes().e_Z) / HEIGHT_RESOLUTION));
      const size_t number_radius_voxels( size_t( UPPER_RADIUS / RADIUS_RESOLUTION));
      const size_t number_angle_voxels( NUMBER_WEDGES);

      // resize data (create tensor of appropriate size, filled with zeros)
      m_Data =
        math::Tensor< double>
       (
         number_height_voxels,
         number_radius_voxels,
         number_angle_voxels,
         double( 0)
       );

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < number_height_voxels; ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( 0); itr_radius < number_radius_voxels; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < number_angle_voxels; ++itr_angle)
          {
            // create cylindrical coordinates
            const coord::CylinderCoordinates cylind_coord
            (
              ( int( itr_length) - int( number_height_voxels / 2)) * m_HeightResolution,
              itr_radius * m_RadiusResolution,
              itr_angle * m_AngleResolution
            );

            // convert to Cartesian coordinates
            linal::Vector3D cart_coord( cylind_coord.GetCartesianCoordinates());

            // transform Cartesian coordinates with Body
            cart_coord.Transform( BODY.GetOrientation());

            // lookup intensity in spline and insert into this density map
            m_Data( itr_length, itr_radius, itr_angle) =
              SPLINE.F( cart_coord.Z(), cart_coord.Y(), cart_coord.X());
          }
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate list of cylindrical density maps from list of bodies (helices) and density map
    //! @param BODIES list of bodies (helices)
    //! @param DENSITY_MAP Cartesian density map
    //! @param HEIGHT_RESOLUTION "voxel size" in height direction
    //! @param RADIUS_RESOLUTION "voxel size" in radius direction
    //! @param NUMBER_WEDGES number of wedges that the density is divided into
    //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
    //! @return return a list of cylindrical density maps
    storage::List< MapCylindrical> MapCylindrical::CalculateCylindricalMaps
    (
      const util::SiPtrList< const assemble::SSEGeometryInterface> &BODIES,
      const Map &DENSITY_MAP,
      const double &HEIGHT_RESOLUTION,
      const double &RADIUS_RESOLUTION,
      const size_t &NUMBER_WEDGES,
      const double &UPPER_RADIUS
    )
    {
      // convert density map into spline
      math::TricubicSpline densitymap_as_spline( DENSITY_MAP.ConvertDensityToSpline());

      // Initialize empty list of cylindrical density maps
      storage::List< MapCylindrical> list_of_cylindrical_density_maps;

      // iterate over all bodies
      for
      (
        util::SiPtrList< const assemble::SSEGeometryInterface>::const_iterator list_itr( BODIES.Begin()),
          list_itr_end( BODIES.End());
        list_itr != list_itr_end;
        ++list_itr
      )
      {
        // use body and spline to construct cylindrical density map
        // insert cylindrical density map into our list of cylindrical density maps
        list_of_cylindrical_density_maps.PushBack
        (
          MapCylindrical
          (
            **list_itr,
            densitymap_as_spline,
            HEIGHT_RESOLUTION,
            RADIUS_RESOLUTION,
            NUMBER_WEDGES,
            UPPER_RADIUS
          )
        );
      }

      // return list of cylindrical density maps
      return list_of_cylindrical_density_maps;
    }

    //! @brief add up all wedges to get a 2D profile (height vs radius - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram x is height, y is radius)
    math::Histogram2D MapCylindrical::TwoDProfileHeightRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all height/radius pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( 0.0, LOWER_RADIUS), // min height and radius
        storage::VectorND< 2, double>( m_HeightResolution, m_RadiusResolution), // bin size for height and radius
        storage::VectorND< 2, size_t>( m_Data.NumberLayers(), ( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution)) // number bins for height and radius
      );
      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // iterate over radius in resolution steps
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // initialize temporary sum of intensities over all angles to 0.0
          double sum_intensities( 0.0);

          // iterate over angles in resolution steps
          // optimize this by using statistics::Sum
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height/radius coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution,
            ( itr_radius + 0.5) * m_RadiusResolution
          );
          // pushback appropriate sum of intensities over the angles into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all radii to get a 2D profile (height vs angle - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram x is height, y is angle)
    math::Histogram2D MapCylindrical::TwoDProfileHeightAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all height/angle pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( 0.0, 0.0), // min height and angle
        storage::VectorND< 2, double>( m_HeightResolution, m_AngleResolution), // binsize for height and angle
        storage::VectorND< 2, size_t>( m_Data.NumberLayers(), m_Data.GetNumberCols()) // number bins for height and angle
      );

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // iterate over angles in resolution steps
        for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
        {
          // initialize temporary sum of intensities over all radii to 0.0
          double sum_intensities( 0.0);

          // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself)
          // optimize this by using statistics::Sum
          for
          (
            size_t itr_radius( lower_radius_itr);
            itr_radius < upper_radius_itr;
            ++itr_radius
          )
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height/angle coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution,
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the angles into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all layers to get a 2D profile (radius vs angle - with intensity coloring)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 2D intensity profile along main axis (in histogram radius is height, y is angle)
    math::Histogram2D MapCylindrical::TwoDProfileRadiusAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram2D of appropriate size to store intensities for all radius/angle pairs
      math::Histogram2D two_d_profile
      (
        storage::VectorND< 2, double>( LOWER_RADIUS, 0.0), // min radius and angle
        storage::VectorND< 2, double>( m_RadiusResolution, m_AngleResolution), // bin size for radius and angle
        storage::VectorND< 2, size_t>( ( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution), m_Data.GetNumberCols()) // number bins for radius and angle
      );
      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over radius in resolution steps
      for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
      {
        // iterate over angle in resolution steps
        for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
        {
          // initialize temporary sum of intensities over all angles to 0.0
          double sum_intensities( 0.0);

          // iterate over height in resolution steps
          // optimize this by using statistics::Sum
          for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct radius/angle coordinates to pushback into 2D profile
          const storage::VectorND< 2, double> values_for_profile
          (
            ( itr_radius + 0.5) * m_RadiusResolution,
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the height into 2D profile
          two_d_profile.PushBack( values_for_profile, sum_intensities);
        }
      }

      //return 2D profile
      return two_d_profile;
    }

    //! @brief add up all wedges and radii to get a 1D profile (height vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileHeight( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram of appropriate size to store intensities for all height values
      const double min_height( 0.0);
      const double binsize_height( m_HeightResolution);
      const size_t number_bins_height( m_Data.NumberLayers());
      math::Histogram one_d_profile( min_height, binsize_height, number_bins_height);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over height in resolution steps
      for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
      {
        // initialize temporary sum of intensities over all radii and angles to 0.0
        double sum_intensities( 0.0);

        // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself)
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_length + 0.5) * m_HeightResolution
          );
          // pushback appropriate sum of intensities over the radii and angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

    //! @brief add up all wedges and layers to get a 1D profile (radius vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {

      BCL_Assert( LOWER_RADIUS > double( 0) && LOWER_RADIUS < UPPER_RADIUS, "minimal radius must be between 0 and upper radius");

      // initialize Histogram of appropriate size to store intensities for all radius values
      const double min_radius( LOWER_RADIUS);
      const double binsize_radius( m_RadiusResolution);
      const size_t number_bins_radius( ( UPPER_RADIUS - LOWER_RADIUS) / m_RadiusResolution);
      math::Histogram one_d_profile( min_radius, binsize_radius, number_bins_radius);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over radius in resolution steps
      for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
      {
        // initialize temporary sum of intensities over all layers and angles to 0.0
        double sum_intensities( 0.0);

        // iterate over height in resolution steps
        for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
        {
          // iterate over angles in resolution steps
          for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct radius coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_radius + 0.5) * m_RadiusResolution
          );
          // pushback appropriate sum of intensities over the height and angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

    //! @brief add up all wedges and layers to get a 1D profile (angle vs intensity)
    //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
    //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
    //! @return 1D intensity profile along main axis
    math::Histogram MapCylindrical::OneDProfileAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const
    {
      // initialize Histogram of appropriate size to store intensities for all angle values
      const double min_angle( 0.0);
      const double binsize_angle( m_AngleResolution);
      const size_t number_bins_angle( m_Data.GetNumberCols());
      math::Histogram one_d_profile( min_angle, binsize_angle, number_bins_angle);

      // initialize minimal and maximal radii for iteration (in terms of Angstroem and bins)
      const size_t lower_radius_itr( size_t( LOWER_RADIUS / m_RadiusResolution));
      const size_t upper_radius_itr( size_t( UPPER_RADIUS / m_RadiusResolution));

      // iterate over angle in resolution steps
      for( size_t itr_angle( 0); itr_angle < m_Data.GetNumberCols(); ++itr_angle)
      {
        // initialize temporary sum of intensities over all radii and height to 0.0
        double sum_intensities( 0.0);

        // iterate over radius in resolution steps (exclude innermost voxels that belong to density rod itself and user-defined outermost voxels))
        for( size_t itr_radius( lower_radius_itr); itr_radius < upper_radius_itr; ++itr_radius)
        {
          // iterate over height in resolution steps
          for( size_t itr_length( 0); itr_length < m_Data.NumberLayers(); ++itr_length)
          {
            sum_intensities += m_Data( itr_length, itr_radius, itr_angle);
          }

          // initialize vector of correct height coordinates to pushback into 1D profile
          const double value_for_profile
          (
            ( itr_angle + 0.5) * m_AngleResolution
          );
          // pushback appropriate sum of intensities over the angles into 1D profile
          one_d_profile.PushBack( value_for_profile, sum_intensities);
        }
      }

      //return 1D profile
      return one_d_profile;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read MapCylindrical from std::istream
    //! @param ISTREAM std::istream from which the density is read
    //! @return std::istream
    std::istream &MapCylindrical::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_HeightResolution, ISTREAM); // height resolution value
      io::Serialize::Read( m_RadiusResolution, ISTREAM); // radius resolution value
      io::Serialize::Read( m_AngleResolution,  ISTREAM); // angle resolution value
      io::Serialize::Read( m_Minimum,          ISTREAM); // Minimal density value
      io::Serialize::Read( m_Maximum,          ISTREAM); // Maximal density value
      io::Serialize::Read( m_Mean,             ISTREAM); // Mean density value
      io::Serialize::Read( m_Body,             ISTREAM); // orientation of density as a body
      io::Serialize::Read( m_Data,             ISTREAM); // density as a tensor

      // end
      return ISTREAM;
    }

    //! @brief write MapCylindrical to std::ostream
    //! @param OSTREAM std::ostream to which the density is written
    //! @param INDENT indent used when writing the density map
    //! @return std::ostream
    std::ostream &MapCylindrical::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HeightResolution,   OSTREAM, INDENT) << '\n'; // height resolution value
      io::Serialize::Write( m_RadiusResolution,   OSTREAM, INDENT) << '\n'; // radius resolution value
      io::Serialize::Write( m_AngleResolution,    OSTREAM, INDENT) << '\n'; // angle resolution value
      io::Serialize::Write( m_Minimum,            OSTREAM, INDENT) << '\n'; // Minimal density value
      io::Serialize::Write( m_Maximum,            OSTREAM, INDENT) << '\n'; // Maximal density value
      io::Serialize::Write( m_Mean,               OSTREAM, INDENT) << '\n'; // Mean density value

      // orientation of density as a body
      io::Serialize::Write( m_Body,      OSTREAM, INDENT) << '\n';

      // density as a tensor does also contain the dimensions of the density map
      io::Serialize::Write( m_Data,      OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
} // namespace bcl
