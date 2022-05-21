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

#ifndef BCL_SCORE_DENSITY_PROFILE_SSE_AGREEMENT_H_
#define BCL_SCORE_DENSITY_PROFILE_SSE_AGREEMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "density/bcl_density.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_simulate_interface.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensityProfileSSEAgreement
    //! @brief this class compares density profiles of density rods to density profiles of placed SSEs
    //! @details given an experimental density map and a set of density rod restraints, 1/2D profiles for those are
    //!          pre calculated. An operator then takes an SSE and a density rod this SSE is placed in, and calculate the
    //!          the profile of this sse within this rod and compares there profiles using CCC
    //!
    //! @see @link example_score_density_profile_sse_agreement.cpp @endlink
    //! @author pereirkn, woetzen
    //! @date Oct 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensityProfileSSEAgreement :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSEGeometryInterface, assemble::SSE, double>
    {

    public:

    ///////////
    // enums //
    ///////////

      //! @enum Profile1DType
      enum Profile1DType
      {
        e_Angle,
        e_Height,
        e_Radius,
        s_NumberProfile1DTypes
      };

      //! @brief Profile1DType as string
      //! @param PROFILE_1D_TYPE the init type
      //! @return the Profile1DType as string
      static const std::string &GetProfile1DTypeString( const Profile1DType &PROFILE_1D_TYPE);

      //! Wrapper enum for Profile1DType
      typedef util::WrapperEnum< Profile1DType, &GetProfile1DTypeString, s_NumberProfile1DTypes> Profile1DTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! map that stores for each body (Geometry) a density Profile (1D)
      storage::Map< util::SiPtr< const assemble::SSEGeometryInterface>, math::Histogram> m_Profiles;

      //! simulator for the density map
      util::ShPtr< density::SimulateInterface> m_Simulator;

      //! height resolution for cylindrical density map calculation
      double m_HeightResolution;

      //! radius resolution for cylindrical density map calculation
      double m_RadiusResolution;

      //! type of profiles used
      Profile1DType m_ProfileType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      // parameters for construction of profile
      static const size_t s_NumberWedges;
      static const double s_UpperRadius;
      static const double s_LowerRadius;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DensityProfileSSEAgreement();

      //! @brief construct from Simulator, density map and body restraints
      //! @param SP_SIMULATOR ShPtr to a simulator
      //! @param DENSITY_MAP the experimental density map
      //! @param BODY_RESTRAINT bodies
      //! @param PROFILE_TYPE 1d profile to be used
      DensityProfileSSEAgreement
      (
        const util::ShPtr< density::SimulateInterface> &SP_SIMULATOR,
        const density::Map &DENSITY_MAP,
        const restraint::Body &BODY_RESTRAINT,
        const Profile1DType PROFILE_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new DensityProfileSSEAgreement
      DensityProfileSSEAgreement *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get scheme
      //! @return scheme
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() takes a BODY and an SSE and compares their profiles
      //! @param BODY the density rod that is occupied by the given SSE
      //! @param SECONDAY_STRUCTURE_ELEMENT the secondary structure that was assigned to that density rod
      //! @return return a double which is the agreement of the rod with the SSE
      double operator()
      (
        const assemble::SSEGeometryInterface &BODY,
        const assemble::SSE &SECONDAY_STRUCTURE_ELEMENT
      ) const;

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

      //! @brief calculate the CCC between two arrays of doubles
      //! @param BEG_PTR_A
      //! @param END_PTR_A
      //! @param BEG_PTR_B
      //! @param END_PTR_B
      //! @return the cross correlation coefficient between range a and b
      static double CCC( const double *BEG_PTR_A, const double *END_PTR_A, const double *BEG_PTR_B, const double *END_PTR_B);

    }; // class DensityProfileSSEAgreement

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_DENSITY_PROFILE_SSE_AGREEMENT_H_
