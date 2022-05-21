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

#ifndef BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_H_
#define BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_neighborhood_interface.h"
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborhoodExposure
    //! @brief a class that scores the exposure of an amino acid in the context of its neighbors
    //!
    //! @see @link example_score_aa_neighborhood_exposure.cpp @endlink
    //! @author woetzen
    //! @date 22.09.2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborhoodExposure :
      public AANeighborhoodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! AAExposure function to be used for calculations
      util::Implementation< assemble::AAExposureInterface> m_AAExposure;

      //! storage map that contains for each membrane environment and aa type a cubic spline
      storage::Map< biol::EnvironmentType, storage::Map< biol::AAType, math::CubicSplineDamped> > m_EnergyFunctionsMembrane;

      //! map for just soluble potentials
      storage::Map< biol::AAType, math::CubicSplineDamped> m_EnergyFunctionsSoluble;

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
      AANeighborhoodExposure();

      //! @brief constructor from ShPtr to aa exposure function
      //! @param SP_AA_EXPOSURE ShPtr to AA exposure function to be used for aa exposure measure and scoring
      AANeighborhoodExposure
      (
        const assemble::AAExposureInterface &SP_AA_EXPOSURE
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AANeighborhoodExposure copied from this one
      AANeighborhoodExposure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
      size_t GetMinimalSequenceSeparation() const;

      //! @brief access to the distance cutoff
      //! @return distance cutoff above which the neighbor does not have influence on the score anymore
      double GetDistanceCutoff() const;

      //! @brief access to the soluble potentials
      //! @return reference to soluble environment potentials
      const storage::Map< biol::AAType, math::CubicSplineDamped> &GetSolublePotentials() const
      {
        return m_EnergyFunctionsSoluble;
      }

      //! @brief access to the membrane potentials for given environment
      //! @param ENVIRONMENT environment for the potentials
      //! @return reference to membrane environment potentials
      const storage::Map< biol::AAType, math::CubicSplineDamped> &GetMembranePotentials
      (
        const biol::EnvironmentType &ENVIRONMENT
      ) const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_name( "AANeighborHoodExposure");
        return s_name;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
      //! @param AA_NEIGHBOR_LIST neighbor list which's center amino acid is scored in the context
      //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
      //! @return the exposure score for this neighbor list
      double operator()
      (
        const assemble::AANeighborList &AA_NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param AA_NEIGHBOR_LIST neighbor list which's center amino acid is scored in the context
      //! @param OSTREAM the std::ostream to be written to
      //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::AANeighborList &AA_NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief read the energy vector for the amino acid neighbor counts in the membrane
      void ReadEnergyVector();

      //! @brief set the members of this object from the given label
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        ReadEnergyVector();
        return true;
      }

    }; // class AANeighborhoodExposure

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_H_
