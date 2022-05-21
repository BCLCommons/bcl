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

#ifndef BCL_SCORE_DATA_SET_PAIRWISE_BIPOLAR_H_
#define BCL_SCORE_DATA_SET_PAIRWISE_BIPOLAR_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetPairwiseBipolar
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_data_set_pairwise_bipolar.cpp @endlink
    //! @author alexanns
    //! @date Jun 29, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetPairwiseBipolar :
      public math::FunctionInterfaceSerializable< restraint::DataSetPairwise, double>
    {

    private:

    //////////
    // data //
    //////////

      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > m_Data;

      //! the scheme of this mutate
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor taking optional scheme
      //! @param SCHEME optional scheme
      DataSetPairwiseBipolar( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking pool
      //! @param SSE_POOL the pool the sse definitions will come from
      //! @param SCHEME optional scheme
      DataSetPairwiseBipolar( const assemble::SSEPool &SSE_POOL, const std::string &SCHEME = GetDefaultScheme());

      //! @brief Clone function
      //! @return pointer to new DataSetPairwiseBipolar
      DataSetPairwiseBipolar *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of a data set
      //! @param DATA data set to be scored
      //! @return the score of the current data set
      double operator()( const restraint::DataSetPairwise &DATA) const;

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

      //! @brief calculates the weight of every residue in helical sses according to its position relative to the ends
      //!        The residue weights go from 0 to 1 along an sse. The end of the sse that has a 0 and the end that has
      //!        a 1 doesn't matter, as long as it is consistent for all sses (which it is). The weights are calculated
      //!        such that the ends of sses that are on the same side of the membrane have the same weight, and the
      //!        weight of residues along sses is the same for each sse as it goes from one side of the membrane to the
      //!        other.
      //! @param SSE_POOL the pool of sses that will be used to calculate the weights
      //! @return map that has for every residue in helical sses the weight associated according to its position relative
      //!         to the end of the sse
      static storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > CalculateResiduePositionWeight( const assemble::SSEPool &SSE_POOL);

      //! @brief calculates the weight of the given residue according to its position away from the end of the given sse
      //! @param SSE the sse whose ends will be used to calculate the weight
      //! @param AA_BASE the residue of interest whose weight will be calculated
      //! @param N_TERMINUS indicates true if nterminus of the sse should have the largest weight - otherwise cterminus
      //! @return double which is the weight of the residue
      static double CalculateWeight( const assemble::SSE &SSE, const biol::AABase &AA_BASE, const bool N_TERMINUS);

    }; // class DataSetPairwiseBipolar

  } // namespace score
} // namespace bcl
#endif // BCL_SCORE_DATA_SET_PAIRWISE_BIPOLAR_H_ 
