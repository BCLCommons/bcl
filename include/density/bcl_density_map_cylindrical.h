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

#ifndef BCL_DENSITY_MAP_CYLINDRICAL_H_
#define BCL_DENSITY_MAP_CYLINDRICAL_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "coord/bcl_coord.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_tensor.h"
#include "math/bcl_math_tricubic_spline.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MapCylindrical
    //! @brief a cylindrical density map.
    //! @details layer is height in cylinder
    //! row is radius in cylinder
    //! col is angle in cylinder
    //!
    //! @see @link example_density_map_cylindrical.cpp @endlink
    //! @author linders, woetzen
    //! @date 08/17/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MapCylindrical :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      assemble::SSEGeometry m_Body;             //!< orientation of density
      double                m_HeightResolution; //!< height resolution of map
      double                m_RadiusResolution; //!< radius resolution of map
      double                m_AngleResolution;  //!< angle resolution of map
      math::Tensor< double> m_Data;             //!< density as a tensor does also contain the dimensions of the density map
      double                m_Minimum;          //!< Minimum intensity
      double                m_Maximum;          //!< Maximum intensity
      double                m_Mean;             //!< Mean    intensity

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      MapCylindrical();

      //! construct MapCylindrical from given Parameters
      //! @param BODY orientation of cylindrical density map as Body
      //! @param DENSITY_MAP the original Euclidean density map from which the cylindrical density map is constructed
      //! @param HEIGHT_RESOLUTION "voxel size" in height direction
      //! @param RADIUS_RESOLUTION "voxel size" in radius direction
      //! @param NUMBER_WEDGES number of wedges that the density is divided into
      //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
      MapCylindrical
      (
        const assemble::SSEGeometryInterface &BODY,
        const Map &DENSITY_MAP,
        const double &HEIGHT_RESOLUTION,
        const double &RADIUS_RESOLUTION,
        const size_t &NUMBER_WEDGES,
        const double &UPPER_RADIUS
      );

      //! construct MapCylindrical from given Parameters
      //! @param BODY
      //! @param SPLINE the spline of the original euclidean density map
      //! @param HEIGHT_RESOLUTION "voxel size" in height direction
      //! @param RADIUS_RESOLUTION "voxel size" in radius direction
      //! @param NUMBER_WEDGES number of wedges that the density is divided into
      //! @param UPPER_RADIUS the maximal radius around main axis that cylindrical density map extents to
      MapCylindrical
      (
        const assemble::SSEGeometryInterface &BODY,
        math::TricubicSpline &SPLINE,
        const double &HEIGHT_RESOLUTION,
        const double &RADIUS_RESOLUTION,
        const size_t &NUMBER_WEDGES,
        const double &UPPER_RADIUS
      );

      //! copy constructor
      MapCylindrical *Clone() const
      {
        return new MapCylindrical( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return number of Elements (i.e. voxels in density map)
      //! @return the size of the tensor as size_t
      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief return Dimensions
      //! @return the dimensions of the tensor (number layers, rows, cols) as storage::VectorND< 3, size_t>
      // important: this is different from cartesian density map implementation!!
      storage::VectorND< 3, size_t> GetDimensions() const
      {
        return storage::VectorND< 3, size_t>( m_Data.NumberLayers(), m_Data.GetNumberRows(), m_Data.GetNumberCols());
      }

      //! @brief return reference to Value (height, radius, angle)
      //! @return the intensity at the tensor value (layer, row, col) as const double
      double const &operator()( const size_t LAYER, const size_t ROW, const size_t COL) const
      {
        return m_Data( LAYER, ROW, COL);
      }

      //! @brief return reference to changeable Value (height, radius, angle)
      //! @return the intensity at the tensor value (layer, row, col) as double
      double &operator()( const size_t LAYER, const size_t ROW, const size_t COL)
      {
        return m_Data( LAYER, ROW, COL);
      }

      //! @brief return reference to changeable Value (height, radius, angle)
      //! @return the intensity at the tensor value (layer, row, col) as double
      double &operator()( const storage::VectorND< 3, size_t> &INDEX)
      {
        return m_Data( INDEX.First(), INDEX.Second(), INDEX.Third());
      }

      //! @brief return reference to Value (height, radius, angle)
      //! @return the intensity at the tensor value (layer, row, col) as const double
      double const &operator()( const storage::VectorND< 3, size_t> &INDEX) const
      {
        return m_Data( INDEX.First(), INDEX.Second(), INDEX.Third());
      }

      //! @brief return minimal value of Map
      //! @return the minimum intensity of the density map as double
      double GetMinimum() const
      {
        return m_Minimum;
      }

      //! @brief return maximal value of Map
      //! @return the maximum intensity of the density map as double
      double GetMaximum() const
      {
        return m_Maximum;
      }

      //! @brief return mean value of Map
      //! @return the average intensity of the density map as double
      double GetMean() const
      {
        return m_Mean;
      }

      //! @brief return const reference to the Tensor density map
      //! @return the density map as math::Tensor< double>
      const math::Tensor< double> &GetData() const
      {
        return m_Data;
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
      static storage::List< MapCylindrical> CalculateCylindricalMaps
      (
        const util::SiPtrList< const assemble::SSEGeometryInterface> &BODIES,
        const Map &DENSITY_MAP,
        const double &HEIGHT_RESOLUTION,
        const double &RADIUS_RESOLUTION,
        const size_t &NUMBER_WEDGES,
        const double &UPPER_RADIUS
      );

      //! @brief add up all wedges to get a 2D profile (height vs radius - with intensity coloring)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 2D intensity profile along main axis (in histogram x is height, y is radius)
      math::Histogram2D TwoDProfileHeightRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

      //! @brief add up all radii to get a 2D profile (height vs angle - with intensity coloring)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 2D intensity profile along main axis (in histogram x is height, y is angle)
      math::Histogram2D TwoDProfileHeightAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

      //! @brief add up all layers to get a 2D profile (radius vs angle - with intensity coloring)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 2D intensity profile along main axis (in histogram radius is height, y is angle)
      math::Histogram2D TwoDProfileRadiusAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

      //! @brief add up all wedges and radii to get a 1D profile (height vs intensity)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 1D intensity profile along main axis
      math::Histogram OneDProfileHeight( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

      //! @brief add up all wedges and layers to get a 1D profile (radius vs intensity)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 1D intensity profile along main axis
      math::Histogram OneDProfileRadius( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

      //! @brief add up all wedges and layers to get a 1D profile (angle vs intensity)
      //! @param LOWER_RADIUS lower radius from which is integrated to calculate the profile
      //! @param UPPER_RADIUS upper radius to which is integrated to calculate the profile
      //! @return 1D intensity profile along main axis
      math::Histogram OneDProfileAngle( const double &LOWER_RADIUS, const double &UPPER_RADIUS) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read MapCylindrical from std::istream
      //! @param ISTREAM std::istream from which the density is read
      //! @return std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write MapCylindrical to std::ostream
      //! @param OSTREAM std::ostream to which the density is written
      //! @param INDENT indent used when writing the density map
      //! @return std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MapCylindrical

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_MAP_CYLINDRICAL_H_ 
