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
#include "chemistry/bcl_chemistry_fragment_connector.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FragmentConnector::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentConnector())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentConnector::FragmentConnector()
    {
    }

    //! @brief Constructor from data members
    //! @param FRAGMENT the base fragment which will be grown by connecting to other fragments
    //! @param VERTICES isomorphism between the base fragment and the molecule which is being assembled
    FragmentConnector::FragmentConnector
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< size_t> &VERTICES
    ) :
      m_Fragment( FRAGMENT),
      m_Vertices( VERTICES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FragmentConnector
    FragmentConnector *FragmentConnector::Clone() const
    {
      return new FragmentConnector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentConnector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief merges the incoming fragment to fragment that is already held by the object over overlapping bonds
    //! @param FRAGMENT the incoming fragment that needs to be merged to the fragment that is already held by this object
    //! @param COMMON_VERTICES atoms that are common between m_Fragment (keys) and incoming fragment (mapped values)
    //! @param ISOMORPHISM isomorphism between incoming fragment and the molecule which is being built
    //! @return true if merge was successful, false otherwise
    bool FragmentConnector::operator()
    (
      const FragmentComplete &FRAGMENT,
      const storage::Map< size_t, size_t> &COMMON_VERTICES,
      const storage::Vector< size_t> &ISOMORPHISM
    )
    {
      storage::Vector< size_t> appended_atoms;
      // merge the provided fragment and the existing fragment that this object holds
      storage::Pair< bool, FragmentComplete> merged_molecule
      (
        MergeFragmentComplete::MergeFragments( m_Fragment, FRAGMENT, COMMON_VERTICES, appended_atoms)
      );

      // if merge was successful then update member data
      if( merged_molecule.First())
      {
        // set assembled fragment to the merged molecule
        m_Fragment = merged_molecule.Second();

        // update the mapping between molecule which is being built and part that has already been built
        for
        (
          storage::Vector< size_t>::const_iterator itr_iso( appended_atoms.Begin()), itr_iso_end( appended_atoms.End());
          itr_iso != itr_iso_end;
          ++itr_iso
        )
        {
          if( !( m_Vertices.Find( ISOMORPHISM( *itr_iso)) < m_Vertices.GetSize()))
          {
            m_Vertices.PushBack( ISOMORPHISM( *itr_iso));
          }
        }
        return true;
      }
      return false;
    }

    //! @brief connects the incoming fragment to fragment that is already held through a bond
    //! @param FRAGMENT the incoming fragment that needs to be added to the fragment that is already held by the object
    //! @param BOND_TYPE the bond that is to be used to connect
    //! @param INDICES_TO_CONNECT atoms of m_Fragment (keys) that need to be connected to atoms of incoming fragment(mapped values)
    //! @param ISOMORPHISM isomorphism between incoming fragment and the molecule which is being built
    //! @return true if merge was successful, false otherwise
    bool FragmentConnector::operator()
    (
      const FragmentComplete &FRAGMENT,
      const ConfigurationalBondType &BOND_TYPE,
      const storage::Pair< size_t, size_t> &INDICES_TO_CONNECT,
      const storage::Vector< size_t> &ISOMORPHISM
    )
    {
      storage::Pair< bool, FragmentComplete> connected_fragment
      (
        MergeFragmentComplete::MergeFragments( m_Fragment, FRAGMENT, BOND_TYPE, INDICES_TO_CONNECT)
      );
      // if connection was successful, then update member variables
      if( connected_fragment.First())
      {
        m_Fragment = connected_fragment.Second();

        for
        (
          storage::Vector< size_t>::const_iterator itr_iso( ISOMORPHISM.Begin()), itr_iso_end( ISOMORPHISM.End());
          itr_iso != itr_iso_end;
          ++itr_iso
        )
        {
          if( !( m_Vertices.Find( *itr_iso) < m_Vertices.GetSize()))
          {
            m_Vertices.PushBack( *itr_iso);
          }
        }
        return true;
      }
      else
      {
        return false;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentConnector::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Fragment, ISTREAM);
      io::Serialize::Read( m_Vertices, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentConnector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Fragment, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Vertices, OSTREAM, INDENT) << '\n';

      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

