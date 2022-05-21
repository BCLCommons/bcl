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

#ifndef BCL_CONTACT_CALCULATE_CORRELATIONS_SM_H_
#define BCL_CONTACT_CALCULATE_CORRELATIONS_SM_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_calculate_correlations_interface.h"
#include "score/bcl_score_aa_assignment_blosum.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

// TODO: Rename!

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CalculateCorrelationsSM
    //! @brief TODO: Calculates a CorrelationMatrix from a MSA which is obtained from the AlignmentInterface class
    //! @details TODO: This class calculates correlated mutations using a similarity matrix which is classically the
    //!                MacLachlan matrix but other matrices can also be used
    //!
    //! @see @link example_contact_calculate_correlations_sm.cpp @endlink
    //! @author teixeipl
    //! @date Jul 31, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CalculateCorrelationsSM :
      public CalculateCorrelationsInterface< biol::AABase>
    {

    private:

    //////////
    // data //
    //////////

        // TODO: Pass in with default make parameter
      static const double s_GapCutoff; // Default portion of sequences with gaps that serves as a cutoff to ignore that column

      //! Similarity matrix to be used for the calculations of correlation
      //! TODO: You should probably make sure that you don't need to do any scaling for your varying types of SMs
      //! TODO: This also should probably be generalized to a type that could take any SM
      static const score::AAAssignmentBLOSUM s_SimilarityMatrix;    // Default is e_BLOSUM_62

      // Map that holds all the weight values for pairs of sequences with keys < key1, key2> being index values for
      // the sequences in their assignment order within the AlignmentInterface and where key1 < key2
      // weights are normalized to sum to 1 and each weight is based on the fraction of non-identical positions
      // between k and l, with larger weights corresponding to less identity
      // making this mutable may not be ideal but it must be seen by GetProbability() which will make the second portion of the code easier to read unless I remove that function...
      mutable storage::Map< storage::Pair< size_t, size_t>, double> m_WeightMap;

    public:

      // single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CalculateCorrelationsSM();

      //! @brief Clone function
      //! @return pointer to new CalculateCorrelationsSM
      CalculateCorrelationsSM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates a correlation matrix if given an AlignmentInterface of type bio::AABase
      //! @param ALIGNMENT_INTERFACE is the MSA representation from which the CM is calculated
      //! @return Returns a correlation matrix of dimensions N by N where N isconst the length of the MSA, containing doubles
      CorrelationMatrix operator()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const;

      //! @brief Checks to see if param CalculateCorrelationsSM object is equal to this one
      //! @param CALCULATE_CORRELATIONS reference to a CalculateCorrelationsSM object
      //! @return bool true for all members are equal and false otherwise
      bool operator==( const CalculateCorrelationsSM &CALCULATE_CORRELATIONS) const;

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

      //! @brief Takes two sequences from the MSA given to CalculateWeights and returns the weight value for that pair
      //! @param SEQUENCE_A from the MSA
      //! @param SEQUENCE_B from the MSA
      //! @return Returns a double weight value
      double GetWeight( const size_t K_SEQ_INDEX, const size_t L_SEQ_INDEX) const;

    }; // class CalculateCorrelationsSM

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_CALCULATE_CORRELATIONS_SM_H_
