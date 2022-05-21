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

#ifndef BCL_SCORE_AA_PAIR_CONTACT_ENERGY_H_
#define BCL_SCORE_AA_PAIR_CONTACT_ENERGY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_pair_distance_interface.h"
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_aa_types.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_3d.h"
#include "math/bcl_math_tricubic_spline.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr.h"
namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairContactEnergy
    //! @brief This is a Function derived template class for scoring AA pair distances
    //!
    //! @see @link example_score_aa_pair_contact_energy.cpp @endlink
    //! @author mendenjl
    //! @date Dec 16, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairContactEnergy :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! scheme to be used in outputting
      std::string m_Scheme;

      typedef storage::Vector< util::ShPtrVector< math::Histogram3D> > AAPairHistogramType;
      typedef storage::VectorND< 3, storage::VectorND< 3, AAPairHistogramType> > SSEPairAAPairHistogramType;

      //! map of energy functions to be used
      util::SiPtr< const SSEPairAAPairHistogramType> m_Histograms;

      //! distance cutoff above which score will always be 0
      double m_DistanceCutoff;

      //! contact probability cutoff
      double m_ContactProbabilityCutoff;

      //! whether to consider loops
      bool m_ConsiderLoops;

      //! whether to consider interface only
      bool m_InterfaceOnly;

      //! interaction and clash instances
      static const util::SiPtr< const util::ObjectInterface> s_InteractionInstance;

    public:

    //////////
    // data //
    //////////

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      AAPairContactEnergy();

      //! @brief virtual copy constructor
      //! @return pointer to a new AAPairContactEnergy object that is copied from this one
      AAPairContactEnergy *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns filename of the histogram being used
      //! @return filename of the histogram being used
      const std::string &GetHistogramFilename() const
      {
        return m_HistogramFileName;
      }

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation in sequence distance
      size_t GetMinimalSequenceSeparation() const
      {
        return 1;
      }

      //! @brief access to the maximal distance that is considered
      //! @return the maximal distance between two amino acids that is scored
      double GetDistanceCutoff() const
      {
        return m_DistanceCutoff;
      }

      //! @brief are amino acids of different chains considered
      //! @return return true, if amino acid pairs of different chains are considered
      bool GetConsiderDifferentChain() const
      {
        return true;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate amino acid pairing potential for given amino acid pair
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param DISTANCE distance between the amino acid pair
      //! @param TYPE_A, TYPE_B SS-type for AMINO_ACID_A/B
      //! @return amino acid pairing potential for given amino acid pair
      double operator()
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const double DISTANCE,
        const biol::SSType &TYPE_A,
        const biol::SSType &TYPE_B
      ) const;

      //! @brief calculate clash score for a protein model
      //! @param MODEL the protein model of interest
      //! @return amino acid pairing potential for given protein
      double operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const biol::SSType &TYPE_A,
        const biol::SSType &TYPE_B,
        std::ostream &OSTREAM
      ) const;

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param MODEL model of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &MODEL,
        std::ostream &OSTREAM
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    private:

      //! @brief read map of amino acid pair energies based on distance from histogram files
      void ReadEnergyFunctionMap();

      //! @brief Get histograms from a particular file. Caches histograms so they need only be read in once
      const SSEPairAAPairHistogramType &GetHistograms( const std::string &FILENAME) const;

    }; // class AAPairContactEnergy

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_PAIR_CONTACT_ENERGY_H_
