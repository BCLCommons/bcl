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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_INVERTER_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_INVERTER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelInverter
    //! @brief Class that all parts of sequence not represented in the given model as SSEs as a inverted protein model
    //! @details This class is developed in conjunction with calculating entropy score terms. For a given protein
    //! model, it first identifies all stretches of AASequence of the chains not represented in the SSEs. Then it looks
    //! at its member m_Sequences which stores extended chain conformation of the sequence for each chain, which is
    //! idealized according similar to a beta-strand, and collects the identified stretches as loop SSEs and returns
    //! them in a new protein model. The returned pseudo-SSEs are cached in member m_SSEs, so that sequence regions
    //! that are usually not represented in the model do not have to be recalculated.
    //!
    //! @see @link example_assemble_protein_model_inverter.cpp @endlink
    //! @author karakam
    //! @date Nov 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelInverter :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! boolean to whether to cache results
      bool m_UseCache;

      //! maximum protein cache size
      size_t m_ModelCacheSize;

      //! maximum coil cache size
      size_t m_CoilCacheSize;

      //! map of extended sequences for given sequences
      mutable storage::Map< util::ShPtr< biol::AASequence>, util::ShPtr< biol::AASequence> > m_ExtendedSequences;

      //! already constructed coils for each given AASequence
      mutable storage::Map< util::ShPtr< biol::AASequence>, util::ShPtrList< SSE> > m_Coils;

      //! already constructed inverted models for each given description key
      mutable storage::List< storage::Pair< std::string, util::ShPtr< ProteinModel> > > m_InvertedModels;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from whether to use cache and a cache size
      //! @param USE_CACHE whether to use cache, by default it is set to false
      //! @param MODEL_CACHE_SIZE Maximum number of models to keep in cache, default is 2
      //! @param COIL_CACHE_SIZE Maximum number of coils to keep in cache, default is 30
      ProteinModelInverter
      (
        const bool USE_CACHE = false,
        const size_t MODEL_CACHE_SIZE = 2,
        const size_t COIL_CACHE_SIZE = 30
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelInverter
      ProteinModelInverter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return whether cache is used
      //! @return whether cache is used
      bool GetUseCache() const
      {
        return m_UseCache;
      }

      //! @brief return model cache size
      //! @return model cache size
      size_t GetModelCacheSize() const
      {
        return m_ModelCacheSize;
      }

      //! @brief return coil cache size
      //! @return coil cache size
      size_t GetCoilCacheSize() const
      {
        return m_CoilCacheSize;
      }

      //! @brief generates the description string for a given chain
      //! @param CHAIN Chain of interest
      //! @return description string
      static std::string GenerateDescription( const Chain &CHAIN);

      //! @brief generates the description string for a given model
      //! @param MODEL ProteinModel of interest
      //! @return description string
      static std::string GenerateDescription( const ProteinModel &MODEL);

      //! @brief returns whether any extended loop with the given begin and end SeqIDs exist for the given AASequence
      //! @param SP_SEQUENCE ShPtr to AASequence of interest
      //! @param BEGIN_SEQ_ID seqid for the first amino acid
      //! @param END_SEQ_ID seqid for the last amino acid
      //! @return whether any extended loop with the given begin and end SeqIDs exist for the given AASequence
      bool DoesContainCoil
      (
        const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
        const size_t BEGIN_SEQ_ID,
        const size_t END_SEQ_ID
      ) const;

      //! @brief reset the cached storage
      void Reset();

    ////////////////
    // operations //
    ////////////////

      //! @brief return inverted model for a given ProteinModel
      //! @param MODEL ProteinModel of interest
      //! @return inverted ProteinModel
      util::ShPtr< ProteinModel> GetInvertedModel( const ProteinModel &MODEL) const;

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

      //! @brief get a cached coil using begin seqid and end seqid
      //! @param SP_SEQUENCE ShPtr to AASequence of the chain inverted
      //! @param BEGIN_SEQ_ID seqid for the first amino acid
      //! @param END_SEQ_ID seqid for the last amino acid
      //! @return coil with specified begin and end seqids
      util::ShPtr< SSE> GetCoilFromCache
      (
        const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
        const size_t BEGIN_SEQ_ID,
        const size_t END_SEQ_ID
      ) const;

      //! @brief returns a chain that contains extended loops containing all amino acids that are not in the given chain
      //! @param CHAIN Chain to be inverted
      //! @return inverted chain
      util::ShPtr< Chain> GetInvertedChainFromCache( const Chain &CHAIN) const;

      //! @brief returns an inverted model
      //! @param MODEL ProteinModel to be inverted
      //! @return ShPtr to inverted model
      util::ShPtr< ProteinModel> GetInvertedModelFromCache( const ProteinModel &MODEL) const;

    public:

      //! @brief static function to generate an extended chain from a given AASequence
      //! @param AA_SEQUENCE AASequence of interest
      //! @return ShPtr to AASequence in the extended chain conformation
      static util::ShPtr< biol::AASequence> ConstructExtendedSequence( const biol::AASequence &AA_SEQUENCE);

      //! @brief static function to generate an inverted chain for a given Chain
      //! @param CHAIN Chain of interest
      //! @return ShPtr to inverted Chain
      static util::ShPtr< Chain> ConstructInvertedChain( const Chain &CHAIN);

      //! @brief static function to generate an inverted model for a given ProteinModel
      //! @param MODEL ProteinMOdel of interest
      //! @return ShPtr to inverted ProteinMOdel
      static util::ShPtr< ProteinModel> ConstructInvertedModel( const ProteinModel &MODEL);

    }; // class ProteinModelInverter

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_MODEL_INVERTER_H_ 
