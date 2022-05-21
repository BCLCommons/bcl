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

#ifndef BCL_QUALITY_RMSD_PREPROCESSOR_H_
#define BCL_QUALITY_RMSD_PREPROCESSOR_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix3x3.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RMSDPreprocessor
    //! @brief RMSDPreprocessor preprocesses a set of coordinates for rapid RMSD calculations with other RMSDPreprocessors
    //!
    //! @see @link example_quality_rmsd_preprocessor.cpp @endlink
    //! @author mendenjl
    //! @date Jun 6, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RMSDPreprocessor :
      public util::ObjectInterface
    {

    private:
    //////////
    // data //
    //////////

      //! square centered norm
      float m_SquareCenteredNorm;

      //! matrix of coordinates recentered, 3xN-coords
      linal::Matrix< float> m_CoordinatesTransposed;

      //! center
      linal::Vector< float> m_Center;

      //! storage space for computing covariance matrix
      mutable linal::Matrix3x3< float> m_Covariance;

      //! bool - whether to recenter (false only if doing real space rmsd)
      bool m_RecenterForSuperimposition;

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
      RMSDPreprocessor() :
        m_SquareCenteredNorm( 0.0),
        m_Center( size_t( 3), float( 0.0)),
        m_RecenterForSuperimposition( false)
      {
      }

      //! @brief constructor from coordinates
      RMSDPreprocessor( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const bool &RECENTER);

      //! @brief virtual copy constructor
      //! @return pointer to new RMSDPreprocessor
      RMSDPreprocessor *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief get the size
      size_t GetSize() const
      {
        return m_CoordinatesTransposed.GetNumberCols();
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT number of indentations
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief determine the transformation matrix to optimally (lowest RMSD_ superimpose two set of coordinates)
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return Transformation matrix that superimposes B onto A
      math::TransformationMatrix3D SuperimposeCoordinates( const RMSDPreprocessor &COORDINATES_B) const;

      //! @brief calculate the rmsd of two sets of coordinates if they are optimally
      //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return the rmsd of the coordinates
      double SuperimposedRMSD( const RMSDPreprocessor &COORDINATES_B) const;

      //! @brief calculate the rmsd of two sets of coordinates if they are optimally
      //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
      //! @param COORDINATES_A set of coordinates A
      //! @param COORDINATES_B set of coordinates B
      //! @return the rmsd of the coordinates
      double RMSD( const RMSDPreprocessor &COORDINATES_B) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
      //! both coordinate sets are translated to the center of mass
      linal::Matrix3x3< float> BuildCovarianceMatrix( const RMSDPreprocessor &COORDINATES_B) const;

    }; // class RMSDPreprocessor

  } // namespace quality
} // namespace bcl

#endif //BCL_QUALITY_RMSD_PREPROCESSOR_H_
