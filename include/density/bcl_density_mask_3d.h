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

#ifndef BCL_DENSITY_MASK_3D_H_
#define BCL_DENSITY_MASK_3D_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_tensor.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mask3d
    //! @brief This class calculates a mask for a area specified by atoms passed to it with the constructor
    //! @details The class calculates the the mask weighing the importance of each point in the grid by regarding its distance
    //! to other atoms, is able to calculate the statistic mean and standard deviation of intensities in an electron
    //! density map, either covered by the mask or in a specified area
    //! the most important function is the calculation of a cross correlation coefficient between two masks with each
    //! point being weighed by the mask
    //!
    //! @see @link example_density_mask_3d.cpp @endlink
    //! @author bitterd
    //! @date Jun 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Mask3d :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      math::Tensor< double>      m_Mask;          //!< actual values for the mask
      storage::VectorND< 3, int> m_Index;         //!< index relative to reference map grid (usually 0,0,0)
      linal::Vector3D             m_Position;      //!< position in the space (should usually be 0,0,0)
      linal::Vector3D             m_GridSpacing;   //!< dimensions of each voxel

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Mask3d();

      //! @brief construct from list of coordinates
      //! @param COORDS a list of coordinates that defines the mask
      //! @param MASKING_DISTANCE distance for sigmoid function to have marginal impact
      //! @param GRID_SPACING length of cube in grid
      //! @param CELL_POSITION position of mask in space
      Mask3d
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDS,
        const double MASKING_DISTANCE,
        const linal::Vector3D &GRID_SPACING,
        const linal::Vector3D &CELL_POSITION
      );

      //! @brief Clone function
      //! @return pointer to new Mask3d
      Mask3d *Clone() const
      {
        return new Mask3d( *this);
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

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the cross correlation between two density maps weighted by mask
      //! @param DENSITY_MAP_EXP experimental density map
      //! @param DENSITY_MAP_SIM simulated density map
      //! @return CCC
      double CrossCorrelationCoefficient
      (
        const Map &DENSITY_MAP_EXP,
        const Map &DENSITY_MAP_SIM
      ) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief calculate mean and standard deviation of over the mask given a density map
      //! is not using the values in the mask, just the box that is defined by the mask
      //! @param DENSITY_MAP density map of interest
      //! @return dataset statistics mean sd
      math::RunningAverageSD< double> CalculateMeanSD( const Map &DENSITY_MAP) const;

      //! @brief calculate mean and standard deviation of over the mask given a density map
      //! is not using the values in the mask, just the box that is defined by the mask
      //! @param DENSITY_MAP
      //! @param OVERLAP
      //! @return dataset statistics mean sd
      math::RunningAverageSD< double> CalculateMeanSDCommonRegion
      (
        const Map &DENSITY_MAP,
        const storage::VectorND< 2, storage::VectorND< 3, size_t> > &OVERLAP
      ) const;

      //! @brief determines the overlap between the mask and two given ranges
      //! @param RANGE_ONE
      //! @param RANGE_TWO
      //! @return returns the common overlap
      static storage::VectorND< 2, storage::VectorND< 3, size_t> > CommonOverlap
      (
        const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_ONE,
        const storage::VectorND< 2, storage::VectorND< 3, size_t> > &RANGE_TWO
      );

      //! @brief reading access to Mask
      //! @return returns Mask
      math::Tensor< double> GetMask() const
      {
        return m_Mask;
      };

      //! @brief reading access to Index
      //! @return returns Index
      storage::VectorND< 3, int> GetIndex() const
      {
        return m_Index;
      }

      //! @brief reading access to grid spacing
      //! @return returns grid spacing
      linal::Vector3D GetGridSpacing() const
      {
        return m_GridSpacing;
      }

      //! @brief reading access to Position
      //! @return returns Position
      linal::Vector3D GetPosition() const
      {
        return m_Position;
      }

      //! @brief determine corners of grid the coordinates fit in
      //! @param COORDS a list of coordinates that defines the grid
      //! @return min coord xyz and max coord xyz
      static storage::VectorND< 2, linal::Vector3D> DetermineGridCorners
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDS
      );

    private:

      //! @brief add margins to corners
      //! @param CORNERS min coord xyz and max coord xyz
      //! @param MARGIN margin to add
      //! @return corners with additional margin
      static storage::VectorND< 2, linal::Vector3D> AddMargin
      (
        const storage::VectorND< 2, linal::Vector3D> &CORNERS,
        const double MARGIN
      );

      //! @brief sigmoid function
      //! @param ARGUMENT_X
      //! @return value of the sigmoide function
      static double Sigmoid( const double &ARGUMENT_X)
      {
        return 1.0 / ( 1.0 + exp( -ARGUMENT_X));
      }

      //! @brief determine the indices ranges that overlap with given density map, relative to this mask
      //! @param DENSITY_MAP overlapping density map
      //! @return pair of start indices and end indices
      storage::VectorND< 2, storage::VectorND< 3, size_t> > OverlappingIndices( const Map &DENSITY_MAP) const;

      //! @brief calculate relative index of given density map to the mask
      //! @param DENSITY_MAP in question
      //! @return relative index
      storage::VectorND< 3, int> RelativeIndex( const Map &DENSITY_MAP) const;

    }; // class Mask3d

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_MASK_3D_H_
