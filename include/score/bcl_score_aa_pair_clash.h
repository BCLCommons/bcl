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

#ifndef BCL_SCORE_AA_PAIR_CLASH_H_
#define BCL_SCORE_AA_PAIR_CLASH_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_pair_distance_interface.h"
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_types.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairClash
    //! @brief This is a scoring function to evaluate whether two amino acids clash
    //! @details The user can specify which types of atoms to use to calculate pairwise distance and the threshold for clash
    //! identification
    //!
    //! @see @link example_score_aa_pair_clash.cpp @endlink
    //! @author karakam
    //! @date 06/15/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairClash :
      public AAPairDistanceInterface,
      public ProteinModel,
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! maximum distance for distance between two amino acids' specified atoms in order to be identified as a clash
      storage::Map< storage::Pair< biol::AAType, biol::AAType>, double>
        m_ShortestObservedDistance;

      //! m_ShortestObservedDistance as a matrix for faster access
      linal::Matrix< double> m_ShortestObservedDistanceMatrix;

      //! width of the sigmoidal function to be used
      double m_SigmoidWidth;

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! distance cutoff above which score will always be 0
      double m_DistanceCutoff;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param SCHEME scheme to be used
      AAPairClash
      (
        const double SIGMOID_WIDTH = 1.0,
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief constructor from a specified histogram file
      //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
      //! @param MINIMAL_SEQUENCE_SEPARATION minimal sequence separation
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param SCHEME scheme to be used
      AAPairClash
      (
        const double SIGMOID_WIDTH,
        const size_t MINIMAL_SEQUENCE_SEPARATION,
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AAPairClash object that is copied from this one
      AAPairClash *Clone() const;

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

      //! @brief create a matrix of all amino acid pair min distances
      //! @return matrix with minimal distances observed for every amino acid pair
      const linal::Matrix< double> &GetShortestObservedDistanceMatrix() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_name( "AAPairClash");
        return s_name;
      }

      //! Get the closest observed distance between two AAs in the PDB
      double GetClosestDistance( const biol::AAType &A, const biol::AAType &B) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate amino acid pairing potential for given amino acid pair
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @return amino acid pairing potential for given amino acid pair
      double operator()( const biol::AABase &AMINO_ACID_A, const biol::AABase &AMINO_ACID_B) const;

      //! @brief calculate amino acid pairing potential for given amino acid pair
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param DISTANCE distance between the amino acid pair
      //! @return amino acid pairing potential for given amino acid pair
      double operator()
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const double DISTANCE
      ) const;

      //! @brief calculate clash score for a protein model
      //! @param MODEL the protein model of interest
      //! @return amino acid pairing potential for given protein
      double operator()( const assemble::ProteinModel &MODEL) const;

      //! @brief calculate clash score for SSEs
      //! @param SSE_A, SSE_B the SSEs to check for clashes
      //! @return clash score
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

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
        std::ostream &OSTREAM
      ) const;

      //! @brief return parameters for data members that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief repulsive term from distance and shortest observed distance
      //! @param DISTANCE the actual distance between amino acids of interest
      //! @param SHORTEST_OBSERVED_DISTANCE shortest distance observed
      double CalculateRepulsiveTerm
      (
        const double DISTANCE,
        const double SHORTEST_OBSERVED_DISTANCE
      ) const;

      //! @brief get the static instance of this class
      static const AAPairClash &GetInstance();

    private:

      //! @brief read map of amino acid pair observed distanced from histogram file
      void ReadDistanceMap();

      //! @brief set the members of this object from the given label
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; //class AAPairClash

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_PAIR_CLASH_H_
