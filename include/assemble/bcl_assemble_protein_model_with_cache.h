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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_WITH_CACHE_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_WITH_CACHE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "descriptor/bcl_descriptor_sequence_interface.h"
#include "signal/bcl_signal_slots.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelWithCache
    //! @brief class that represents a biological protein structure as a model within the bcl.
    //! @details An adaptor for a protein model to a descriptor::SequenceInterface
    //!
    //! @see @link example_assemble_protein_model_with_cache.cpp @endlink
    //! @author mendenjl
    //! @date Jan 07, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelWithCache :
      public ProteinModel,
      public descriptor::SequenceInterface< biol::AABase>,
      public signal::Slots
    {

    protected:

    //////////
    // data //
    //////////

      //! this is the complete sequence belonging to this model
      util::SiPtrVector< biol::AABase> m_AAs;

      //! whether to only retrieve AAs that have defined coordinates
      bool m_RequireCoordinates;

      //! pointer to the fragment complete used for this class
      mutable util::ShPtr< chemistry::AAFragmentComplete> m_FragmentComplete;

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
      ProteinModelWithCache();

      //! @brief construct from a protein model
      //! @param MODEL protein model of interest
      //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
      ProteinModelWithCache( const ProteinModel &MODEL, const bool &REQUIRE_COORDINATES);

      //! @brief copy constructor
      //! @param ORIGINAL model with cache to copy
      ProteinModelWithCache( const ProteinModelWithCache &ORIGINAL);

      //! @brief copy constructor
      //! @return pointer a new ProteinModelWithCache copied from this model
      virtual ProteinModelWithCache *Clone() const;

      //! @brief hard copy constructor
      //! @return a ProteinModelWithCache with chains hard copied from that model
      virtual ProteinModelWithCache *HardCopy() const;

      //! @brief empty copy constructor
      //! @return a ProteinModelWithCache that is empty
      virtual ProteinModelWithCache *Empty() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns whether coordinates are required
      bool GetRequiresCoordinates() const
      {
        return m_RequireCoordinates;
      }

      //! @brief return the length of the sequence in question
      //! @return the length of the sequence in question
      size_t GetSize() const;

      //! @brief get the iterator for the sequence
      //! @return the iterator for the sequence
      iterate::Generic< const biol::AABase> GetIterator() const;

      //! @brief get a non-constant iterator for the sequence
      //! @return the non-constant iterator for the sequence
      iterate::Generic< biol::AABase> GetIteratorNonConst();

      //! @brief get the protein model, represented as a molecule
      //! @note this will only be constructed if when this function is called
      const chemistry::AAFragmentComplete &GetChemicalRepresentation() const;

      //! @brief Reset the cache
      virtual void ResetCache() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param PROTEIN_MODEL ProteinModelWithCache to be copied
      //! @return This model after all members are assigned to values from PROTEIN_MODEL
      virtual ProteinModelWithCache &operator =( const ProteinModelWithCache &PROTEIN_MODEL);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ProteinModelWithCache from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProteinModelWithCache to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief update amino acid pointers based on the new model
      void UpdateAAPtrs( const ProteinModel &MODEL);

    }; // class ProteinModelWithCache

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PROTEIN_MODEL_WITH_CACHE_H_
