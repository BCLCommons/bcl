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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_HINGE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_HINGE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_move_interface.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEPairHinge
    //! @brief applies provided moves to domains collected around hinge point if provided
    //! @details This Mutate Class uses its collector m_Collector to collect a SiPtrList of SSEs and converts them to a domain
    //! and applies the provided move by m_Move to the newly created domain. If m_Locator is provided, this move will
    //! be relative to the hinge point located by m_Locator, otherwise the hinge point will be defined by the domain
    //! itself
    //!
    //! @see @link example_fold_mutate_protein_model_sse_pair_hinge.cpp @endlink
    //! @author karakam, woetzen
    //! @date May 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateProteinModelSSEPairHinge :
      public math::MutateInterface< assemble::ProteinModel>
    {
    private:

    //////////
    // data //
    //////////

      //! collector for pairs of sses
      util::Implementation< find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface> > m_Collector;

      //! moves the sse
      util::ShPtr< coord::MoveInterface> m_Move;

      //! boolean to decide whether to move the hinge or not
      bool m_MoveHinge;

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
      MutateProteinModelSSEPairHinge();

      //! @brief constructor from a CollectorInterface, MoveInterface and a scheme
      //! @param COLLECTOR function that chooses the sses
      //! @param MOVE function that performs the move on the sse
      //! @param MOVE_HINGE whether the hinge should be also moved
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEPairHinge
      (
        const find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface> &COLLECTOR,
        const coord::MoveInterface &MOVE,
        const bool MOVE_HINGE = true,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEPairHinge>()
      );

      //! @brief constructor from ShPtrs to CollectorInterface, MoveInterface and a scheme
      //! @param SP_COLLECTOR ShPtr to function that chooses the sses
      //! @param SP_MOVE ShPtr to function that performs the move on the sse
      //! @param MOVE_HINGE whether the hinge should be also moved
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEPairHinge
      (
        const util::ShPtr< find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface> > &SP_COLLECTOR,
        const util::ShPtr< coord::MoveInterface> &SP_MOVE,
        const bool MOVE_HINGE = true,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEPairHinge>()
      );

      //! @brief clone
      MutateProteinModelSSEPairHinge *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "MutateProteinModelSSEPairHinge");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Applies transformations to domains around hinge points.");
        serializer.AddInitializer
        (
          "collector",
          "collector for domains",
          io::Serialization::GetAgent( &m_Collector)
        );
        // serializer.AddInitializer
        // (
        //   "move",
        //   "move object for SSEs",
        //   io::Serialization::GetAgent( &m_Move)
        // );
        serializer.AddInitializer
        (
          "move_hinge",
          "move hinges",
          io::Serialization::GetAgent( &m_MoveHinge)
        );
        return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL ProteinModel which will be mutated
      //! @return MutateResult ProteinModel after mutation
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

    }; // class MutateProteinModelSSEPairHinge

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_PAIR_HINGE_H_
