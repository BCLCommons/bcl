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

#ifndef BCL_ASSEMBLE_PICK_PROTEIN_MODEL_CONFORMATION_RANDOM_H_
#define BCL_ASSEMBLE_PICK_PROTEIN_MODEL_CONFORMATION_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_pick_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickProteinModelConformationRandom
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_pick_protein_model_conformation_random.cpp @endlink
    //! @author alexanns
    //! @date Apr 15, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickProteinModelConformationRandom :
      public find::PickInterface< util::SiPtr< const ProteinModel>, util::SiPtrList< const ProteinModel> >
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickProteinModelConformationRandom();

      //! @brief Clone function
      //! @return pointer to new PickProteinModelConformationRandom
      PickProteinModelConformationRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Pick returns a random SSE from the domain argument
      //! @param ENSEMBLE model which the locator will pick from
      //! @return returns SiPtr to protein model conformation from protein model conformation ensemble
      util::SiPtr< const ProteinModel> Pick( const util::SiPtrList< const ProteinModel> &ENSEMBLE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class PickProteinModelConformationRandom

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PICK_PROTEIN_MODEL_CONFORMATION_RANDOM_H_ 
