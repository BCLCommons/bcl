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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_protein_model_multiplier.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_transformer.h"
#include "math/bcl_math_function_cached.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelMultiplier::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelMultiplier())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelMultiplier::ProteinModelMultiplier() :
      m_ChainMultipliers(),
      m_NumberMultimers( util::GetUndefinedSize_t()),
      m_Orientation()
    {
    }

    //! @brief construct from triplets of chain ids and transformations
    //! @param TRANSFORMATIONS vector of current chain id, new chain id and transformation
    //! @param PROTEIN_MODEL original protein model
    //! @param CACHE bool whether to cache the SSE transformer
    ProteinModelMultiplier::ProteinModelMultiplier
    (
      const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > &TRANSFORMATIONS,
      const ProteinModel &PROTEIN_MODEL,
      const bool CACHE
    ) :
      m_ChainMultipliers(),
      m_NumberMultimers( util::GetUndefinedSize_t()),
      m_Orientation()
    {
      BuildChainMultipliers( TRANSFORMATIONS, PROTEIN_MODEL, CACHE);
    }

    //! @brief construct from an axis of symmetry and number of subunits
    //! @param SYMMETRY_AXIS axis of symmetry
    //! @param SUBUNITS number of subunits
    //! @param PROTEIN_MODEL original protein model
    //! @param CACHE bool whether to cache the SSE transformer
    ProteinModelMultiplier::ProteinModelMultiplier
    (
      const linal::Vector3D &SYMMETRY_AXIS,
      const size_t SUBUNITS,
      const ProteinModel &PROTEIN_MODEL,
      const bool CACHE
    ) :
      m_ChainMultipliers(),
      m_NumberMultimers( SUBUNITS),
      m_Orientation()
    {
      BuildChainMultipliers( GetChainData( SYMMETRY_AXIS, SUBUNITS, PROTEIN_MODEL.GetChainIDs()), PROTEIN_MODEL, CACHE);
    }

    //! @brief construct from an axis of symmetry and number of subunits
    //! @param SYMMETRY_AXIS axis of symmetry
    //! @param SUBUNITS number of subunits
    //! @param PROTEIN_MODEL original protein model
    //! @param DIHEDRAL_AXIS optional secondary rotation axis for dihedral symmetry
    //! @param CACHE bool whether to cache the SSE transformer
    ProteinModelMultiplier::ProteinModelMultiplier
    (
      const linal::Vector3D &SYMMETRY_AXIS,
      const linal::Vector3D &DIHEDRAL_AXIS,
      const size_t SUBUNITS,
      const ProteinModel &PROTEIN_MODEL,
      const bool CACHE
    ) :
      m_ChainMultipliers(),
      m_NumberMultimers( SUBUNITS),
      m_Orientation()
    {
      BuildChainMultipliers
      (
        GetChainData( SYMMETRY_AXIS, SUBUNITS, PROTEIN_MODEL.GetChainIDs(), DIHEDRAL_AXIS),
        PROTEIN_MODEL,
        CACHE
      );
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelMultiplier
    ProteinModelMultiplier *ProteinModelMultiplier::Clone() const
    {
      return new ProteinModelMultiplier( *this);
    }

    //! @brief make a copy of this multiplier that also copies the ShPtrs in the set of chain multipliers
    //! @return new copy of this
    util::ShPtr< ProteinModelMultiplier> ProteinModelMultiplier::HardCopy() const
    {
      storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan> data_copy;
      // iterate through the set to make copies
      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          itr( m_ChainMultipliers.Begin()), itr_end( m_ChainMultipliers.End());
        itr != itr_end;
        ++itr
      )
      {
        data_copy.Insert( ( *itr)->HardCopy());
      }
      util::ShPtr< ProteinModelMultiplier> copy( new ProteinModelMultiplier());
      copy->m_ChainMultipliers = data_copy;
      copy->m_NumberMultimers = m_NumberMultimers;
      copy->m_Orientation = m_Orientation;

      return copy;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelMultiplier::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the orientation of the object
    //! @return orientation
    linal::Vector3D ProteinModelMultiplier::GetAxis( const coord::Axis &AXIS) const
    {
      return m_Orientation.GetAxis( AXIS);
    }

    //! @brief gets the transformation matrices for this multiplier
    //! @return the transformation matrices for this multiplier
    storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >
    ProteinModelMultiplier::GetTransformationMatrices() const
    {
      // initialize matrices
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > matrices;

      // iterate through the chain multipliers
      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          multiplier_itr( m_ChainMultipliers.Begin()),
          multiplier_itr_end( m_ChainMultipliers.End());
        multiplier_itr != multiplier_itr_end; ++multiplier_itr
      )
      {
        matrices.PushBack
        (
          storage::Triplet< char, char, math::TransformationMatrix3D>
          (
            ( *multiplier_itr)->GetInitialChainID(),
            ( *multiplier_itr)->GetNewChainID(),
            *( ( *multiplier_itr)->GetTransformationMatrix())
          )
        );
      }

      // end
      return matrices;
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D ProteinModelMultiplier::GetCenter() const
    {
      return m_Orientation.GetOrigin();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief transform the multimer
    //! @param TRANSFORMATION transformation to be applied
    void ProteinModelMultiplier::Transform( const math::TransformationMatrix3D &TRANSFORMATION)
    {
      // iterate through the chain multipliers
      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          multiplier_itr( m_ChainMultipliers.Begin()),
          multiplier_itr_end( m_ChainMultipliers.End());
        multiplier_itr != multiplier_itr_end; ++multiplier_itr
      )
      {
        // get the current transformation
        util::ShPtr< math::TransformationMatrix3D> sp_chain_transformation
        (
          ( *multiplier_itr)->GetTransformationMatrix()
        );

        // move back to the original location, apply the chain transformation, then apply the passed transformation
        math::TransformationMatrix3D transform( math::Inverse( TRANSFORMATION));
        transform( *sp_chain_transformation);
        transform( TRANSFORMATION);
        *sp_chain_transformation = transform;
      }

      m_Orientation( TRANSFORMATION);
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void ProteinModelMultiplier::Translate( const linal::Vector3D &TRANSLATION)
    {
      Transform( math::TransformationMatrix3D( TRANSLATION));
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void ProteinModelMultiplier::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      Transform( math::TransformationMatrix3D( ROTATION_MATRIX_3D));
    }

    //! @brief mapping of chain id mapping from origianl to the new chain id
    //! @return table containing the matrix number, the original chainid, the chainid it will have in the model
    storage::Table< char> ProteinModelMultiplier::ChainIDMapping() const
    {
      storage::Table< char> mapping( storage::Vector< std::string>::Create( "orig", "new"));
      size_t count( 0);

      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          itr( m_ChainMultipliers.Begin()), itr_end( m_ChainMultipliers.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        mapping.InsertRow
        (
          util::Format()( count),
          storage::Vector< char>::Create( ( *itr)->GetInitialChainID(), ( *itr)->GetNewChainID())
        );
      }

      return mapping;
    }

    //! @brief gets a map of original chain id to target chain ids (as a string)
    //! @return map of original chain id to target chain ids (as a string)
    storage::Map< char, std::string> ProteinModelMultiplier::GetTargetChains() const
    {
      // initialize map
      storage::Map< char, std::string> chains;

      // iterate over chain multipliers
      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          itr( m_ChainMultipliers.Begin()), itr_end( m_ChainMultipliers.End());
        itr != itr_end; ++itr
      )
      {
        // add chain id char to string
        chains[ ( *itr)->GetInitialChainID()].push_back( ( *itr)->GetNewChainID());
      }

      // end
      return chains;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns protein model after transforming all chains
    //! @param PROTEIN_MODEL protein model to be replicated
    //! @return protein model after transforming all chains
    ProteinModel ProteinModelMultiplier::operator ()( const ProteinModel &PROTEIN_MODEL) const
    {
      // initialize vector of chains
      util::ShPtrVector< Chain> chains;

      // iterate through the chain multipliers
      for
      (
        storage::Set< util::ShPtr< ChainMultiplier>, ChainMultiplierLessThan>::const_iterator
          multiplier_itr( m_ChainMultipliers.Begin()), multiplier_itr_end( m_ChainMultipliers.End());
        multiplier_itr != multiplier_itr_end; ++multiplier_itr
      )
      {
        // get the new chain
        util::ShPtr< Chain> new_chain
        (
          ( *multiplier_itr)->operator ()( *( PROTEIN_MODEL.GetChain( ( *multiplier_itr)->GetInitialChainID())))
        );

        // push back into the vector
        chains.PushBack( new_chain);
      }

      // sort the chains by chain id
      chains.Sort( ChainLessThan());
      ProteinModel new_model( chains);

      // set the protein model data
      util::ShPtr< ProteinModelData> pmd( PROTEIN_MODEL.GetProteinModelData()->HardCopy());
      pmd->Replace( ProteinModelData::e_Multiplier, util::ShPtr< ProteinModelMultiplier>());
      new_model.SetProteinModelData( pmd);

      // end
      return new_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelMultiplier::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChainMultipliers, ISTREAM);
      io::Serialize::Read( m_NumberMultimers, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelMultiplier::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainMultipliers, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberMultimers, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief wrap an SSE transformer into a FunctionCached object
    //! @param SSE_TRANSFORMER SSE transformer to be wrapped
    //! @param CACHE bool whether to cache the SSE transformer
    //! @return FunctionCached object
    util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > ProteinModelMultiplier::WrapCacheSSETransformer
    (
      const util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > > &SSE_TRANSFORMER,
      const bool CACHE
    )
    {
      // does function need to be wrapped
      if( !CACHE)
      {
        // return the function
        return SSE_TRANSFORMER;
      }
      else
      {
        // wrap function into cache function
        util::ShPtr< math::FunctionCached< SSE, util::ShPtr< SSE> > > sp_cache_function
        (
          new math::FunctionCached< SSE, util::ShPtr< SSE> >
          (
            SSE_TRANSFORMER,
            &SSE::GetDestructorSignal
          )
        );

        // add signal handler for coordinate changes
        sp_cache_function->AddSignalHandlerForArgument( &SSE::GetCoordinateChangeSignal);

        // end
        return sp_cache_function;
      }
    }

    //! @brief converts axis, subunit, and chain id information into triplets used by the constructor
    //! @param SYMMETRY_AXIS axis of symmetry
    //! @param SUBUNITS number of subunits
    //! @param CHAIN_IDS chain ids present in original model
    //! @param DIHEDRAL_AXIS secondary rotation axis for dihedral symmetry
    //! @return transformation information for each chain multiplier
    storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > ProteinModelMultiplier::GetChainData
    (
      const linal::Vector3D &SYMMETRY_AXIS,
      const size_t SUBUNITS,
      const std::string &CHAIN_IDS,
      const linal::Vector3D &DIHEDRAL_AXIS
    )
    {
      // make sure the number of subunits is positive
      BCL_Assert
      (
        SUBUNITS != 0 && SUBUNITS != util::GetUndefined< size_t>(),
        "Invalid number of symmetric subunits for protein model multiplier"
      );

      // initialize chain data
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > chain_data;

      // make sure at least one chain id is given
      if( CHAIN_IDS.empty())
      {
        return chain_data;
      }

      // identity transformation
      const math::TransformationMatrix3D transform_identity;

      // iterate through the chains in the original protein model
      for
      (
        std::string::const_iterator chain_itr( CHAIN_IDS.begin()), chain_itr_end( CHAIN_IDS.end());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // pushback the identity transformation
        chain_data.PushBack
        (
          storage::Triplet< char, char, math::TransformationMatrix3D>( *chain_itr, *chain_itr, transform_identity)
        );
      }

      // get the new chain id
      char new_chain_id( CHAIN_IDS[ CHAIN_IDS.length() - 1]);

      // increment through the number of subunits
      for( size_t i( 1); i != SUBUNITS; ++i)
      {
        // calculate the transformation
        const math::TransformationMatrix3D transform
        (
          math::RotationMatrix3D( SYMMETRY_AXIS, double( i) * 2.0 * math::g_Pi / double( SUBUNITS))
        );

        // iterate through the original chains
        for
        (
          std::string::const_iterator chain_itr( CHAIN_IDS.begin()), chain_itr_end( CHAIN_IDS.end());
          chain_itr != chain_itr_end; ++chain_itr
        )
        {
          // increment to next available chain id
          ++new_chain_id;

          // pushback the transformation
          chain_data.PushBack
          (
            storage::Triplet< char, char, math::TransformationMatrix3D>( *chain_itr, new_chain_id, transform)
          );
        }
      }

      // if this is dihedral symmetry
      if( DIHEDRAL_AXIS != linal::Vector3D())
      {
        // copy the chain data
        storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > orig_chain_data( chain_data);

        // initialize subunit counter
        size_t subunit( 0);

        // iterate through the transformations
        for
        (
          storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >::const_iterator
            itr( orig_chain_data.Begin()), itr_end( orig_chain_data.End());
          itr != itr_end; ++itr, ++subunit
        )
        {
          // calculate the transformation
          math::TransformationMatrix3D transform( itr->Third());
          transform( math::RotationMatrix3D( DIHEDRAL_AXIS, math::g_Pi));

          // iterate through the original chains
          for
          (
            std::string::const_iterator chain_itr( CHAIN_IDS.begin()), chain_itr_end( CHAIN_IDS.end());
            chain_itr != chain_itr_end; ++chain_itr
          )
          {
            // increment to next available chain id
            ++new_chain_id;

            // pushback the transformation
            chain_data.PushBack
            (
              storage::Triplet< char, char, math::TransformationMatrix3D>( *chain_itr, new_chain_id, transform)
            );
          }
        }
      }

      // end
      return chain_data;
    }

    //! @brief builds the member data from the transformations and protein model
    //! @param TRANSFORMATIONS vector of current chain id, new chain id and transformation
    //! @param PROTEIN_MODEL original protein model
    //! @param CACHE bool whether to cache the SSE transformer
    void ProteinModelMultiplier::BuildChainMultipliers
    (
      const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > &TRANSFORMATIONS,
      const ProteinModel &PROTEIN_MODEL,
      const bool CACHE
    )
    {
      // initialize coords
      storage::Vector< linal::Vector3D> coords;

      // iterate through the vector
      for
      (
        storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >::const_iterator
          transformation_itr( TRANSFORMATIONS.Begin()), transformation_itr_end( TRANSFORMATIONS.End());
        transformation_itr != transformation_itr_end; ++transformation_itr
      )
      {
        // get the chain
        const util::ShPtr< Chain> &sp_chain( PROTEIN_MODEL.GetChain( transformation_itr->First()));
        if( !sp_chain.IsDefined())
        {
          BCL_MessageCrt
          (
            "no such chain in model: " + std::string( 1, transformation_itr->First())
          );
          continue;
        }
        // construct the new sequence for the chain at the position of the chain
        util::ShPtr< biol::AASequence> sp_sequence
        (
          sp_chain->GetSequence()->HardCopy()
        );
        sp_sequence->SetChainID( transformation_itr->Second());
        sp_sequence->Transform( transformation_itr->Third());

        // get the transformation matrix
        const util::ShPtr< math::TransformationMatrix3D> sp_transformation_matrix( transformation_itr->Third().Clone());

        // add to coords
        coords.PushBack( transformation_itr->Third().GetOrigin());

        // create a chain multiplier and add it to the set
        m_ChainMultipliers.Insert
        (
          util::ShPtr< ChainMultiplier>
          (
            new ChainMultiplier
            (
              WrapCacheSSETransformer
              (
                util::ShPtr< util::FunctionInterface< SSE, util::ShPtr< SSE> > >
                (
                  new SSETransformer( sp_sequence, sp_transformation_matrix)
                ),
                CACHE
              ),
              transformation_itr->First(),
              sp_transformation_matrix,
              sp_sequence
            )
          )
        );
      }

      // set the orientation
      m_Orientation = math::TransformationMatrix3D( coord::CenterOfMass( util::SiPtrVector< const linal::Vector3D>( coords.Begin(), coords.End())));
    }

  } // namespace assemble
} // namespace bcl
