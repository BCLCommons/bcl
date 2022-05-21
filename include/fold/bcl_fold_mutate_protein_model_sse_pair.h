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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEPair
    //! @brief transforms an SSE relative to a second SSE
    //! @details Using a Collector, pairs of SSEs within a protein model are collected. A randomly a SSE pair is selected.
    //! In that pair one sse selected to be mutated, the other one gives the relative position. Along the shortest
    //! connection between those two bodies a translation and rotation is applied.
    //!
    //! @see @link example_fold_mutate_protein_model_sse_pair.cpp @endlink
    //! @author woetzen, karakam
    //! @date May 24, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSEPair :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! collector for pairs of sses
      util::ShPtr
      <
        find::CollectorInterface
        <
          storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >,
          assemble::DomainInterface
        >
      > m_Collector;

      //! min translation alonbg shortest connection
      math::Range< double> m_TranslationRange;

      //! min and max rotation around shortest connection
      math::Range< double> m_RotationRange;

      //! scheme
      std::string m_Scheme;

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
      MutateProteinModelSSEPair();

      //! @brief construct from Collector, max translation and rotation
      //! @param SP_COLLECTOR ShPtr to collector of SSE pairs
      //! @param MAX_TRANSLATION maximum translation allowed
      //! @param MAX_ROTATION maximum rotation allowed
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEPair
      (
        const util::ShPtr
        <
          find::CollectorInterface
          <
            storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
          >
        > &SP_COLLECTOR,
        const double MAX_TRANSLATION,
        const double MAX_ROTATION,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEPair>()
      );

      //! @brief construct from Collector, min and max rotation and translation
      //! @param SP_COLLECTOR ShPtr to collector of SSE pairs
      //! @param TRANSLATION_RANGE range that specifies min and max translations
      //! @param ROTATION_RANGE range that specifies min and max rotations
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEPair
      (
        const util::ShPtr
        <
          find::CollectorInterface
          <
            storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
          >
        > &SP_COLLECTOR,
        const math::Range< double> &TRANSLATION_RANGE,
        const math::Range< double> &ROTATION_RANGE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEPair>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSEPair
      MutateProteinModelSSEPair *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // MutateProteinModelSSEPair

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_H_ 
