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

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_pick_atom_by_element.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickAtomByElement::s_Instance
    (
      util::Enumerated< find::PickInterface< util::SiPtr< const AtomConformationalInterface>, util::SiPtrList< const AtomConformationalInterface> > >::AddInstance( new PickAtomByElement())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    //! @details will pick from carbons as the default element type
    PickAtomByElement::PickAtomByElement() :
      m_ElementType( GetElementTypes().e_Carbon)
    {
    }

    //! @brief constructor from element type
    //! @param ELEMENT_TYPE element which should be picked
    PickAtomByElement::PickAtomByElement( const ElementType &ELEMENT_TYPE) :
      m_ElementType( ELEMENT_TYPE)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    PickAtomByElement *PickAtomByElement::Clone() const
    {
      return new PickAtomByElement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickAtomByElement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PickAtomByElement::GetAlias() const
    {
      static const std::string s_alias( "PickAtomByElement");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PickAtomByElement::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Picks a random atom of a specified element from a set.");
      serializer.AddInitializer
      (
        "element",
        "type of the element to pick",
        io::Serialization::GetAgent( &m_ElementType)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Picks a random atom of the appropriate type
    //! @param ATOMS is the SiPtrList which provides the pool of Atom to pick from
    //! @return returns SiPtr to a random Atom in the list
    util::SiPtr< const AtomConformationalInterface>
      PickAtomByElement::Pick( const util::SiPtrList< const AtomConformationalInterface> &ATOMS) const
    {
      BCL_Assert( !ATOMS.IsEmpty(), "atom list is empty!!!")

      // atom_subset is the subset of ATOMS which match the desired parameters (i.e. match the element type)
      util::SiPtrList< const AtomConformationalInterface> atom_subset;
      for
      (
        util::SiPtrList< const AtomConformationalInterface>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetElementType() == m_ElementType)
        {
          atom_subset.PushBack( *itr);
        }
      }

      BCL_Assert( !atom_subset.IsEmpty(), "no atoms of appropriate type found!!!")

      // return
      return *random::GetGlobalRandom().Iterator( atom_subset.Begin(), atom_subset.End(), atom_subset.GetSize());
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
