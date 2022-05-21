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
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "linal/bcl_linal_vector_operations.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from an input stream
    //! @param ISTREAM input stream
    //! @param ADD_H whether to add hydrogens to the ensemble automatically if the format is not e_BCL
    //! @param RANGE the range of small molecules to load from the stream
    FragmentEnsemble::FragmentEnsemble
    (
      std::istream &ISTREAM,
      const sdf::HydrogenHandlingPref &H_PREF,
      const math::Range< size_t> &RANGE,
      const sdf::NeutralizationPref &NEUTRALIZATION
    )
    {
      ReadMoreFromMdl( ISTREAM, H_PREF, RANGE, NEUTRALIZATION);
    }

    //! @brief Clone function
    //! @return pointer to new FragmentEnsemble
    FragmentEnsemble *FragmentEnsemble::Clone() const
    {
      return new FragmentEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read additional molecules into the ensemble
    //! @param ISTREAM input stream, reads in SDF format
    //! @param ADD_H whether to add hydrogens to the ensemble automatically if the format is e_BCL
    //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
    void FragmentEnsemble::ReadMoreFromMdl
    (
      std::istream &ISTREAM,
      const sdf::HydrogenHandlingPref &H_PREF,
      const math::Range< size_t> &RANGE,
      const sdf::NeutralizationPref &NEUTRALIZATION
    )
    {
      // get the original size
      const size_t original_size( GetSize());
      // close the borders on the range for ease of use
      math::Range< size_t> closed_range( RANGE.CloseBorders());

      // intended # to read; max is necessary in case the range is the range of size_t
      const size_t n_to_read( std::max( closed_range.GetWidth(), closed_range.GetWidth() + size_t( 1)));

      // read in from the stream until we reach the end of the file or the last index
      for( FragmentFeed feed( ISTREAM, H_PREF, n_to_read, closed_range.GetMin(), NEUTRALIZATION); feed.NotAtEnd(); ++feed)
      {
        m_FragmentEnsemble.PushBack( *feed);
      }

      BCL_MessageVrb( "finished reading ensemble with " + util::Format()( GetSize() - original_size) + " molecules.");
    }

    //! @brief shuffle the ensemble
    void FragmentEnsemble::Shuffle()
    {
      std::vector< storage::List< FragmentComplete>::iterator> ens_list;
      ens_list.reserve( m_FragmentEnsemble.GetSize());
      for( auto itr( m_FragmentEnsemble.Begin()), itr_end( m_FragmentEnsemble.End()); itr != itr_end; ++itr)
      {
        ens_list.push_back( itr);
      }
      std::random_shuffle( ens_list.begin(), ens_list.end(), random::GetRandomSizeT);
      storage::List< FragmentComplete> new_frag_ensemble;
      for( auto itr( ens_list.begin()), itr_end( ens_list.end()); itr != itr_end; ++itr)
      {
        new_frag_ensemble.InternalData().splice( new_frag_ensemble.End(), m_FragmentEnsemble.InternalData(), *itr);
      }
      std::swap( new_frag_ensemble, m_FragmentEnsemble);
    }

    //! @brief Remove all hydrogens
    void FragmentEnsemble::RemoveH()
    {
      for
      (
        storage::List< FragmentComplete>::iterator
          itr_mols( m_FragmentEnsemble.Begin()), itr_mols_end( m_FragmentEnsemble.End());
        itr_mols != itr_mols_end;
        ++itr_mols
      )
      {
        itr_mols->RemoveH();
      }
    }

    //! @brief Add hydrogens
    void FragmentEnsemble::SaturateWithH()
    {
      for
      (
        storage::List< FragmentComplete>::iterator
          itr_mols( m_FragmentEnsemble.Begin()), itr_mols_end( m_FragmentEnsemble.End());
        itr_mols != itr_mols_end;
        ++itr_mols
      )
      {
        itr_mols->SaturateWithH();
      }
    }

    //! @brief sort ensemble based on increasing value property
    //! @param PROPERTY property of interest which is stored on the molecule
    //! @return a vector with molecules sorted on the basis of property in increasing order
    void FragmentEnsemble::Sort
    (
      const std::string PROPERTY_NAME
    )
    {
      std::vector< std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> > ens_list;
      ens_list.reserve( m_FragmentEnsemble.GetSize());
      for( auto itr( m_FragmentEnsemble.Begin()), itr_end( m_FragmentEnsemble.End()); itr != itr_end; ++itr)
      {
        ens_list.push_back( std::make_pair( itr->GetStoredProperties().GetMDLPropertyAsVector( PROPERTY_NAME), itr));
      }
      std::stable_sort
      (
        ens_list.begin(),
        ens_list.end(),
        [](
            const std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> &A,
            const std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> &B
        )->bool { return A.first < B.first;}
      );
      storage::List< FragmentComplete> new_frag_ensemble;
      for( auto itr( ens_list.begin()), itr_end( ens_list.end()); itr != itr_end; ++itr)
      {
        new_frag_ensemble.InternalData().splice( new_frag_ensemble.End(), m_FragmentEnsemble.InternalData(), itr->second);
      }
      std::swap( new_frag_ensemble, m_FragmentEnsemble);
    }

    //! @brief sort ensemble based on increasing value property
    //! @param PROPERTY property of interest which is stored on the molecule
    //! @return a vector with molecules sorted on the basis of property in increasing order
    void FragmentEnsemble::Sort( const storage::Vector< std::string> PROPERTY_NAMES)
    {
      std::vector< std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> > ens_list;
      ens_list.reserve( m_FragmentEnsemble.GetSize());
      size_t res_size( 0);
      for( auto itr( m_FragmentEnsemble.Begin()), itr_end( m_FragmentEnsemble.End()); itr != itr_end; ++itr)
      {
        storage::Vector< float> props;
        props.AllocateMemory( res_size);
        for( auto itr_p( PROPERTY_NAMES.Begin()), itr_p_end( PROPERTY_NAMES.End()); itr_p != itr_p_end; ++itr_p)
        {
          auto vec( itr->GetStoredProperties().GetMDLPropertyAsVector( *itr_p));
          props.Append( storage::Vector< float>( vec.Begin(), vec.End()));
        }
        res_size = std::max( res_size, props.GetSize());
        ens_list.push_back( std::make_pair( linal::Vector< float>( props), itr));
      }
      std::stable_sort
      (
        ens_list.begin(),
        ens_list.end(),
        [](
            const std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> &A,
            const std::pair< linal::Vector< float>, storage::List< FragmentComplete>::iterator> &B
        )->bool { return A.first < B.first;}
      );
      storage::List< FragmentComplete> new_frag_ensemble;
      for( auto itr( ens_list.begin()), itr_end( ens_list.end()); itr != itr_end; ++itr)
      {
        new_frag_ensemble.InternalData().splice( new_frag_ensemble.End(), m_FragmentEnsemble.InternalData(), itr->second);
      }
      std::swap( new_frag_ensemble, m_FragmentEnsemble);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &FragmentEnsemble::WriteMDL( std::ostream &OSTREAM) const
    {
      // iterate through all molecules in the given ensemble
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_mols( m_FragmentEnsemble.Begin()),
          itr_mols_end( m_FragmentEnsemble.End());
        itr_mols != itr_mols_end;
        ++itr_mols
      )
      {
        itr_mols->WriteMDL( OSTREAM);
      }
      return OSTREAM;
    }

    std::istream &FragmentEnsemble::Read( std::istream &ISTREAM)
    {
      ReadMoreFromMdl( ISTREAM, sdf::e_Maintain);
      return ISTREAM;
    }

    std::ostream &FragmentEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      WriteMDL( OSTREAM);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
