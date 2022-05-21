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

#ifndef BCL_FOLD_MUTATE_DOMAIN_SSE_PAIR_TRIM_H_
#define BCL_FOLD_MUTATE_DOMAIN_SSE_PAIR_TRIM_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "math/bcl_math_mutate_interface.h"
#include "pdb/bcl_pdb_factory.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDomainSSEPairTrim
    //! @brief Trims a number of residues two sses such that after trimming the
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_fold_mutate_domain_sse_pair_trim.cpp @endlink
    //! @author alexanns
    //! @date Sep 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDomainSSEPairTrim :
      public math::MutateInterface< assemble::Domain>
    {

    private:

    //////////
    // data //
    //////////

      //! the scheme for this mutate
      std::string m_Scheme;

      //! minimum SSE sizes to be allowed when shrinking
      storage::Map< biol::SSType, size_t> m_MinSSESizes;

      //! the total number of residues that will possibly be trimmed
      size_t m_NumberResiduesToTrim;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateDomainSSEPairTrim();

      //! @brief constructor taking member variable parameters
      //! @param SCHEME the scheme for this mutate
      //! @param MIN_SSE_SIZES minimum SSE sizes to be allowed when shrinking
      //! @param NUMBER_RESIS_TO_TRIM the total number of residues that will possibly be trimmed
      explicit MutateDomainSSEPairTrim
      (
        const size_t NUMBER_RESIS_TO_TRIM,
        const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES = pdb::Factory::GetCommandlineSSETypeMinSizes(),
        const std::string &SCHEME = GetStaticClassName< MutateDomainSSEPairTrim>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDomainSSEPairTrim
      MutateDomainSSEPairTrim *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a domain
      //! @param DOMAIN domain which will be mutated
      //! @return MutateResult with the mutated domain
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &THIS_DOMAIN) const;

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

      storage::VectorND< 2, util::ShPtr< assemble::SSE> > GetBestTrimmedSSEs
      (
        const assemble::SSE &SSE_N_TERMINAL, const assemble::SSE &SSE_C_TERMINAL
      ) const;

      util::ShPtr< assemble::SSE> TrimSSE
      (
        const assemble::SSE &SSE, const biol::AASequenceFlexibility::SequenceDirection &SEQUENCE_DIRECTION,
        const size_t NUM_RESI_TO_REMOVE
      ) const;

      double CalculateDistancePerResidue
      (
        const assemble::SSE &SSE_N_TERMINAL, const assemble::SSE &SSE_C_TERMINAL
      ) const;

    }; // class MutateDomainSSEPairTrim

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_DOMAIN_SSE_PAIR_TRIM_H_ 
