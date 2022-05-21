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

#ifndef BCL_DENSITY_PROTEIN_AGREEMENT_LIKELIHOOD_H_
#define BCL_DENSITY_PROTEIN_AGREEMENT_LIKELIHOOD_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_density_map.h"
#include "bcl_density_protein_agreement_interface.h"
#include "bcl_density_simulate_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_log_likelihood.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAgreementLikelihood
    //! @brief calculates "density agreement score" between a protein model and a density map
    //! @details Based on the log likelihood of the local cross correlations between residues and masked map, this class
    //! calculates a "density agreement" measure that judges how well a protein model agrees with a given density map
    //!
    //! @see @link example_density_protein_agreement_likelihood.cpp @endlink
    //! @author bitterd
    //! @date Jun 29, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAgreementLikelihood :
      public ProteinAgreementInterface
    {

    private:

    //////////
    // data //
    //////////

      //! scores the likelihood of seeing a given cross correlation coefficient
      math::LogLikelihood m_LogLikelihood;

      //! use all atoms of residue plus neighboring residues, otherwise just CA
      bool m_HighResolution;

      //! set of atomtypes to consider
      storage::Set< biol::AtomType> m_AtomTypes;

      //! density map to calculate the CCC against
      util::SiPtr< const Map> m_DensityMap;

      //! density map simulator to be used
      util::ShPtr< SimulateInterface> m_Simulator;

      //! scheme of density agreement
      std::string m_Scheme;

      //! @brief low res masking distance for Mask3D using CA only
      static const double s_LowResolutionMaskingDistance;

      //! @brief high res masking distance for Mask3D using all side chain atoms
      static const double s_HighResolutionMaskingDistance;

      //! @brief fitted function, for mean CCC by chance
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > m_MeanFit;

      //! @brief fitted function, for sd CCC by chance
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > m_SdFit;

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
      ProteinAgreementLikelihood();

      //! @brief constructor from atom types
      //! @param HIGH_RESOLUTION use all atoms of residue plus neighboring residues, otherwise just CA
      //! @param ATOM_TYPES atom types to consider
      //! @param MEAN_CCC_RESOLUTION_FUNCTION
      //! @param SD_CCC_RESOLUTION_FUNCTION
      ProteinAgreementLikelihood
      (
        const bool HIGH_RESOLUTION,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        const math::FunctionInterfaceSerializable< double, double> &MEAN_CCC_RESOLUTION_FUNCTION,
        const math::FunctionInterfaceSerializable< double, double> &SD_CCC_RESOLUTION_FUNCTION
      );

      //! @brief Clone function
      //! @return pointer to new ProteinAgreementLikelihood
      ProteinAgreementLikelihood *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to the density simulator
      //! @return the density simulator used
      const util::ShPtr< SimulateInterface> &GetSimulator() const;

      //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
      //! @param SP_SIMULATOR ShPtr to SimulatInterface
      void SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR);

      //! @brief access to the density used for agreement calculation
      //! @return SiPtr to the density
      const util::SiPtr< const Map> &GetDensity() const;

      //! @brief set the density used for agreement calculation
      //! @param SP_DENSITY SiPtr to the density map
      void SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the standard deviation between the protein model and the given density
      //! @param PROTEIN_MODEL protein of interest
      //! @return score for the density likelihood
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

      //! @brief calculate low resolution score, only considering CA atoms
      //! @param PROTEIN_MODEL protein of interest
      //! @return score for the density likelihood
      double LowResolutionScore( const assemble::ProteinModel &PROTEINMODEL) const;

      //! @brief calculate high resolution score, considering all atoms of each amino acid sidechain plus the two
      //! neighboring residues
      //! @param PROTEIN_MODEL protein of interest
      //! @return score for the density likelihood
      double HighResolutionScore( const assemble::ProteinModel &PROTEINMODEL) const;

    private:

      //! @brief string for the default scheme
      const std::string &GetDefaultScheme() const;

      //! @brief set scheme
      //! @details set scheme from the default scheme and atom types
      void SetScheme();

    }; // class ProteinAgreementLikelihood

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_PROTEIN_AGREEMENT_LIKELIHOOD_H_ 
