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

#ifndef BCL_CHEMISTRY_SUBSTITUENT_CONFORMATIONAL_H_
#define BCL_CHEMISTRY_SUBSTITUENT_CONFORMATIONAL_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SubstituentConformational
    //! @brief holds the SubstituentConformationals of an atom.
    //! @brief SubstituentConformationals can be ordered according to Cahn-Ingold-Prelog priorty rules using std::less
    //! The SubstituentConformational performs breadth-first searches as necessary when compared with other SubstituentConformationals
    //! Likewise, when determining e.g. stereocenters, the SubstituentConformational will remain of minimal size unless the molecule
    //! has symmetry such that it requires more detail to determine priority.
    //! Known limitation:
    //!  Emergent stereocenters are not considered by operator <
    //!    Emergent stereocenters are formed when two otherwise identical SubstituentConformationals differ only in the chirality of
    //!    their SubstituentConformational atoms
    //!
    //! @see @link example_chemistry_substituent_conformational.cpp @endlink
    //! @author mendenjl
    //! @date Jan 26, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API SubstituentConformational :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Sorted vectors of atoms at increasing distances
      //! mutable to allow operator < to continue the breadth 1st search if necessary to determine priority
      mutable storage::List< util::SiPtrVector< const AtomConformationalInterface> > m_Substituents;

      //! This vector keeps track of atoms that were ghosted in the last round
      //! Ghost atoms arise in CHIP rules from double, triple, and aromatic bonds
      //! @see http://en.wikipedia.org/wiki/Cahn%E2%80%93Ingold%E2%80%93Prelog_priority_rules
      mutable util::SiPtrVector< const AtomConformationalInterface>                  m_PriorGhosts;

      //! Atoms whose bonds have been expanded already
      mutable storage::Set< util::SiPtr< const AtomConformationalInterface> >        m_ExpandedAtoms;

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
      SubstituentConformational();

      //! @brief construct a SubstituentConformational from an ATOM and it's MOLECULE
      SubstituentConformational( const AtomConformationalInterface &ATOM);

      //! @brief virtual copy constructor
      SubstituentConformational *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief get the atom at the root of this search
      //! @return the atom at the root of this search
      const util::SiPtr< const AtomConformationalInterface> &GetRootAtom() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief test if one SubstituentConformational is less than another
      //! @param RHS right hand SubstituentConformationals
      //! @return true if this has higher priority than RHS
      bool operator <( const SubstituentConformational &RHS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT nr of indentations
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief tell whether or not the SubstituentConformational is complete, e.g. whether it has any leaves that can be expanded
      //! @return true if there are no further atoms to consider
      bool IsComplete() const;

      //! @brief add the SubstituentConformationals at the next level of depth
      //! @return the SubstituentConformationals at the next level of description
      const util::SiPtrVector< const AtomConformationalInterface> &Expand() const;

    }; // class SubstituentConformational

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SUBSTITUENT_CONFORMATIONAL_H_

