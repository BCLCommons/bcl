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

#ifndef BCL_SCORE_PROTEIN_MODEL_MEMBRANE_TOPOLOGY_H_
#define BCL_SCORE_PROTEIN_MODEL_MEMBRANE_TOPOLOGY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelMembraneTopology
    //! @brief Penalizes topologies that would force loop/coil regions to pass through the membrane
    //! @details Uses an SSE pool that provides the expected transmembrane helical SSEs. These need to be unbroken,
    //!          continuous transmembrane helices. The pool must be non-overlapping. These criteria are asserted in the
    //!          constructor. The expected TM sses are implicitely considered to be numbered such that every other
    //!          sse has it's n-terminus on the same side of the membrane. The sses from the model are matched with
    //!          expected TM sses. The n-termini of model sses are checked to make sure they are on the same side of
    //!          the membrane (i.e. sign of z-coordinate is the same), for ones that should be on the same side
    //!          according to the expected TM sse they correspond to. The score is incremented for each pair that do
    //!          not have n-termini on the same side, although they should.
    //!
    //! @see @link example_score_protein_model_membrane_topology.cpp @endlink
    //! @author alexanns
    //! @date Jun 14, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelMembraneTopology :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! locator to find the transmembrane domain of the protein model
      util::ShPtr< find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> > m_TransmembraneDomain;

      //! the pool of expected SSEs making up the transmembrane domain
      assemble::SSEPool m_SSEPool;

      //! the scheme for this score
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief gives the flag that allows the input of the filename of the pool containing the expectd tm helices
      //! @return shptr to flag with parameter for providing the expected transmembrane helices
      static const util::ShPtr< command::FlagInterface> &GetFlagExpectedTransmembraneHelicesPoolFile();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelMembraneTopology();

      //! @brief constructor taking member variables
      //! @param LOCATOR the locator that will be used to locate the transmembrane domains of the protein model
      //! @param POOL the pool of sses identifying the transmembrane sses
      //! @param SCHEME the scheme of this score
      ProteinModelMembraneTopology
      (
        const util::ShPtr< find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> > &LOCATOR,
        const assemble::SSEPool &POOL,
        const std::string &SCHEME = GetStaticClassName< ProteinModelMembraneTopology>()
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelMembraneTopology
      ProteinModelMembraneTopology *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const;

      //! @brief get score type
      //! @return score type
      ProteinModel::Type GetType() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize this score from the flag
      //! @return true if initialization was successful - e.g. the flag was set
      bool InitializeFromFlag();

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an PROTEIN_MODEL and returning a t_ResultType object
      //! @param PROTEIN_MODEL Protein Model to be used to evaluate the function
      //! @return function value of the given argument
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

    private:

      //! @brief associate the expected tm sses coming from the sse pool with the sses from the located tm domain
      //! @param TRANSMEMBRANE_DOMAIN
      //! @param EXPECTED_TM_SEGMENTS
      //! @return map with the expected tm sse as the key and the set of associated model tm sses as the value.
      //!         there can be multiple model sses associated with an expected tm sse since in the model the whole tm
      //!         segment might be broken
      static storage::Map
      <
        util::SiPtr< const assemble::SSE>,
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
        assemble::SSELessThan
      > AssociateModelTMSSEsWithExpectedTMSegments
      (
        const assemble::Domain &TRANSMEMBRANE_DOMAIN, const assemble::SSEPool &EXPECTED_TM_SEGMENTS
      );

      //! @brief For a side of the membrane, gathers the coordinates of the n-termini of SSEs that should be on that side
      //!        The side is implicitly defined by the fact that all expected TM sses are provided, and every other SSE
      //!        should have its n-terminus on the same side
      //! @param EXPECTED_AND_MODEL_TM_SSES the assignment of an expected TM sse with corresponding sses from the model
      //! @param ALL_EXPECTED_TM_SSES the list of all expected SSEs
      //! @return map which has a size_t representing a side of the membrane, and the corresponding coordinates from
      //!         model sses assigned to expected tm sses whose n-termini should be on the same side of the membrane
      //!         The char is the chain id of the sse, because if the sses aren't in the same chain they should not be
      //!         compared together
      static storage::Map< size_t, storage::Vector< storage::Pair< double, char> > >
      GetExpectedMembraneSideAndModelSSETerminusCoordinates
      (
        const storage::Map
        <
          util::SiPtr< const assemble::SSE>,
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>,
          assemble::SSELessThan
        > &EXPECTED_AND_MODEL_TM_SSES,
        const util::SiPtrVector< const assemble::SSE> &ALL_EXPECTED_TM_SSES,
        const bool N_TERMINUS
      );

      //! @brief scores the agreement of the actual arrangement of sses with the expected arrangement
      //! @param TOPOLOGY map that implicitly holds the expected and model topologies. The size_t indicates a side of
      //!        the membrane, the vector of doubles is the coordinates of the nterminus of sses that should be on that
      //!        side of the membrane.
      //! @return double which is the score of the agreement of the model sse topology with expected topology
      static double ScoreTopology
      (
        const storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > &TOPOLGY
      );

      //! @brief scores the agreement of the actual arrangement of sses with the expected arrangement
      //!        Penalizes occurances of ntermini being on the same side of the membrane although they shouldn't be,
      //!        since their arrangement would force a loop to go through the membrane
      //! @param TOPOLOGY map that implicitly holds the expected and model topologies. The size_t indicates a side of
      //!        the membrane, the vector of doubles is the coordinates of the nterminus of sses that should be on that
      //!        side of the membrane.
      //! @return double which is the score of the agreement of the model sse topology with expected topology
      static double ScoreTopologyOpposingSides
      (
        const storage::Map< size_t, storage::Vector< storage::Pair< double, char> > > &TOPOLGY
      );

    }; // class ProteinModelMembraneTopology

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_MEMBRANE_TOPOLOGY_H_
