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

#ifndef BCL_CHEMISTRY_SUB_FRAGMENT_H_
#define BCL_CHEMISTRY_SUB_FRAGMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_const_graph.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SubFragment
    //! @brief a helper class for the FragmentMolecule class
    //! @details stores a
    //!
    //! @see @link example_chemistry_sub_fragment.cpp @endlink
    //! @author kothiwsk
    //! @date May 03, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SubFragment :
      public util::ObjectInterface
    {

    private:

      //! Vector for storing vertex isomorphism between this sub fragment and the original molecule that was passed to
      //! FragmentMolecule class
      storage::Vector< size_t> m_ThisToNode;

      //! Vector for storing vertex isomorphism between this sub fragment and its immediate parent molecule which is a
      //! fragment of the original molecule which was passed
      storage::Vector< size_t> m_ThisToParent;

      //! A set representation of m_ThisToNode
      storage::Set< size_t> m_ThisToNodeSet;

      //! A set representation of m_ThisToParent
      storage::Set< size_t> m_ThisToParentSet;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      SubFragment();

      //! @brief constructor given atoms with conformation and molecule configuration
      //! @params MOLECULE  molecule from which sub fragment, which is the molecule itself, has to be created
      SubFragment
      (
        const FragmentComplete &MOLECULE
      );

      //! @brief constructor given atoms with conformation and molecule configuration
      //! @params PARENT_FRAGMENT SubFragment from which this sub_fragment has to be created
      //! @params SUB_INDICES indices of atoms of SubFragment from which this sub fragment has to be created
      SubFragment
      (
        const SubFragment &PARENT_FRAGMENT,
        const storage::Vector< size_t> &SUB_INDICES
      );

      //! @brief Clone function
      //! @return pointer to new SubFragment
      SubFragment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns vector of indices from this subfragment to the original molecule
      //! @return the vector of indices from this subfragment to the original molecule
      const storage::Vector< size_t> &GetThisToNode() const;

      //! @brief returns vector of indices from this subfragment to the parent molecule
      //! @return the vector of indices from this subfragment to the parent molecule
      const storage::Vector< size_t> &GetThisToParent() const;

      //! @brief returns set of indices from this subfragment to the original molecule
      //! @return the set of indices from this subfragment to the original molecule
      const storage::Set< size_t> &GetThisToNodeSet() const;

      //! @brief returns set of indices from this subfragment to the parent molecule
      //! @return the set of indices from this subfragment to the parent molecule
      const storage::Set< size_t> &GetThisToParentSet() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

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
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class SubFragment

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SUB_FRAGMENT_H_
