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

#ifndef BCL_CHEMISTRY_PICK_ATOM_BY_ELEMENT_H_
#define BCL_CHEMISTRY_PICK_ATOM_BY_ELEMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_element_types.h"
#include "find/bcl_find_pick_interface.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickAtomByElement
    //! @brief Class is used for picking a random atom from a subset based on element
    //!
    //! @author morettr
    //! @date 03/27/2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickAtomByElement :
      public find::PickInterface
      <
        util::SiPtr< const AtomConformationalInterface>,
        util::SiPtrList< const AtomConformationalInterface>
      >
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! element type which will be picked
      ElementType                                            m_ElementType;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      PickAtomByElement();

      //! @brief constructor from element type
      //! @param ELEMENT_TYPE element which should be picked
      PickAtomByElement( const ElementType &ELEMENT_TYPE);

      //! virtual copy constructor
      PickAtomByElement *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return element type for picker
      //! @return the element being picked
      const ElementType &GetElementType() const
      {
        return m_ElementType;
      }

      void SetElementType( const ElementType &ELEMENT_TYPE)
      {
        m_ElementType = ELEMENT_TYPE;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Picks a random atom of the appropriate type
      //! @param ATOMS is the SiPtrList which provides the pool of Atom to pick from
      //! @return returns SiPtr to a random Atom in the list of the appropriate type
      util::SiPtr< const AtomConformationalInterface>
        Pick( const util::SiPtrList< const AtomConformationalInterface> &ATOMS) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PickAtomByElement

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_PICK_ATOM_BY_ELEMENT_H_
