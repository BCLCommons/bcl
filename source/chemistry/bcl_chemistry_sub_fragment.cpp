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
#include "chemistry/bcl_chemistry_sub_fragment.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SubFragment::s_Instance
    (
      GetObjectInstances().AddInstance( new SubFragment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    SubFragment::SubFragment()
    {
    }

    //! @brief constructor given atoms with conformation and molecule configuration
    //! @params MOLECULE  molecule from which sub fragment, which is the molecule itself, has to be created
    SubFragment::SubFragment
    (
      const FragmentComplete &MOLECULE
    )
    {
      for( size_t i( 0); i < MOLECULE.GetNumberAtoms(); ++i)
      {
        m_ThisToNode.PushBack( i);
        m_ThisToParent.PushBack( i);
      }
      m_ThisToNodeSet = m_ThisToParentSet = storage::Set< size_t>( m_ThisToNode.Begin(), m_ThisToNode.End());
    }

    //! @brief constructor given atoms with conformation and molecule configuration
    //! @params SUB_FRAGMENT SUB_FRAGMENT from which the sub_fragment has to be created
    //! @params SUB_INDICES indices of SUB_FRAGMENT from which this sub fragment has to be created
    SubFragment::SubFragment
    (
      const SubFragment &PARENT_FRAGMENT,
      const storage::Vector< size_t> &SUB_INDICES
    ) :
      m_ThisToParent( SUB_INDICES),
      m_ThisToParentSet( SUB_INDICES.Begin(), SUB_INDICES.End())
    {
      const storage::Vector< size_t> &parent_node( PARENT_FRAGMENT.GetThisToNode());
      m_ThisToNode.AllocateMemory( SUB_INDICES.GetSize());

      for( size_t i( 0), sz( SUB_INDICES.GetSize()); i < sz; ++i)
      {
        m_ThisToNode.PushBack( parent_node( SUB_INDICES( i)));
      }
      m_ThisToNodeSet.InsertElements( m_ThisToNode.Begin(), m_ThisToNode.End());
    }

    //! @brief Clone function
    //! @return pointer to new SubFragment
    SubFragment *SubFragment::Clone() const
    {
      return new SubFragment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SubFragment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns map between indices of Node to Parent
    //! @return the map between indices of Node to Parent
    const storage::Vector< size_t> &SubFragment::GetThisToNode() const
    {
      return m_ThisToNode;
    }

    //! @brief returns map between indices of Parent to This
    //! @return the map between indices of Parent to This
    const storage::Vector< size_t> &SubFragment::GetThisToParent() const
    {
      return m_ThisToParent;
    }

    //! @brief returns map between indices of Node to Parent
    //! @return the map between indices of Node to Parent
    const storage::Set< size_t> &SubFragment::GetThisToNodeSet() const
    {
      return m_ThisToNodeSet;
    }

    //! @brief returns map between indices of Parent to This
    //! @return the map between indices of Parent to This
    const storage::Set< size_t> &SubFragment::GetThisToParentSet() const
    {
      return m_ThisToParentSet;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SubFragment::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SubFragment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
