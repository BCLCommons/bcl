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

#ifndef BCL_RESTRAINT_SAS_DEBYE_H_
#define BCL_RESTRAINT_SAS_DEBYE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_debye_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasDebye
    //! @brief Uses the Debye formula to calculate I from the form factors for a protein model
    //! @details I(q) = (sum i=1 to M)*(sum j=1 to M) Fi(q)*Fj(q) * sin(q*rij)/(q*rij)
    //!
    //! @see @link example_restraint_saxs_debye.cpp @endlink
    //! @author putnamdk
    //! @date May 5, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasDebye :
      public SasDebyeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! bool value to represent loops that are not present in the protein model
      bool m_ShouldApproximateLoops;

      //! whether to determine the norm factor with regula falsi (true) or pythagorean approximation (false)
      bool m_DetermineAnalyticNormFactor;

      //! Excluded Volume Parameter -  // m_BaseAtomType->GetElementType()->GetStructureFactor() + H * Hform_factor - solv
      //return m_Vacuo->operator ()( Qvalue) - m_Mass * std::exp( math::Sqr( Qvalue) * m_SurfaceArea); user provided input to adjust excluded volume in the form factor calculation
      double m_ExcludedVolumeParameter;

      //! Hydration Shell Parameter - user provided input to adjust the hydration shell in the form factor calculation
      double m_HydrationShellParameter;

      //! bool value to control side chain approximation functionality
      bool m_ShouldApproximateSideChains;

      //! bool value to use Sans implementation of Debye formula if false, use SAXS
      bool m_UseSans;

      //! value between 0 and 1 to control the degree of deuteration of the solvent
      double m_DeuteriumExchangeParameter;

      //! Storage Vector of all Atoms
      mutable storage::Vector< std::string> m_CompleteAllAtoms;
      mutable storage::Vector< linal::Vector3D> m_AtomicCoordinates;
      mutable storage::Vector< double> m_AtomSasa;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Constructor from members
      //! @param LOOPS bool value to represent loops that are not present in the protein model
      //! @param USE_REGULA_FALSI_APPROXIMATION whether to determine the norm factor with regula falsi (true) or
      //!        pythagorean approximation (false)
      //! @param EXCLUDED_VOLUME_PARAMETER value to tune excluded volume in form factor calculation
      //! @param HYDRATION_SHELL_PARAMETER value to tune hydration shell in form factor calculation
      SasDebye
      (
        const bool LOOPS = false,
        const bool USE_REGULA_FALSI_APPROXIMATION = false,
        double EXCLUDED_VOLUME_PARAMETER = 1.0,
        double HYDRATION_SHELL_PARAMETER = 0.0,
        const bool SIDE_CHAIN_APPROXIMATION = false,
        const bool USE_SANS = false,
        double DEUTERIUM_EXCHANGE_PARAMETER = 0.0,
        util::ShPtr< storage::Vector< SasScatteringPoint> > REDUCED_EXP_DATA =
        util::ShPtr< storage::Vector< SasScatteringPoint> >()
      );

      //! @brief Clone function
      //! @return pointer to new SasDebye
      SasDebye *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      const storage::Vector< std::string> &GetCompleteAllAtoms() const
      {
        return m_CompleteAllAtoms;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Get the positions and form factor functions for all the atoms in the protein model
      //! @param MODEL the protein model of interest
      //! @return a vector containing pairs of position and form factor function
      storage::Vector
      <
        storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
      > GetAtomsAndFormFactors( const assemble::ProteinModel &PROTEIN_MODEL) const;

      storage::Vector< std::string> &GetAtomGroups()
      {
        return m_CompleteAllAtoms;
      }

      storage::Vector< linal::Vector3D> &GetCoordinates()
      {
        return m_AtomicCoordinates;
      }

      storage::Vector< double> &GetSASAPoint()
      {
        return m_AtomSasa;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief overloaded () operator to calculate Intensity from Q
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns Intensity for given Q for both the experimental and calculated data
      SasExperimentalAndCalculatedData operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class SasDebye

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_DEBYE_H_
