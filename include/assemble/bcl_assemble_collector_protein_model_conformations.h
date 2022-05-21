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

#ifndef BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATIONS_H_
#define BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATIONS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorProteinModelConformations
    //! @brief collects all conformations of a protein model into a list
    //! @details All conformations representing a protein model are gathered into a list. The current model of focus
    //!          is optionally collected depending on the member variable.
    //!
    //! @see @link example_assemble_collector_protein_model_conformations.cpp @endlink
    //! @author alexanns
    //! @date May 3, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorProteinModelConformations :
      public find::CollectorInterface< util::SiPtrList< const ProteinModel>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! if true the current conformation held by the protein model will be considered with the other conformations
      bool m_ConsiderCurrentConformation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorProteinModelConformations();

      //! @brief constructor taking parameters
      //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
      CollectorProteinModelConformations( const bool CONSIDER_CURRENT);

      //! @brief Clone function
      //! @return pointer to new CollectorProteinModelConformations
      CollectorProteinModelConformations *Clone() const;

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

      //! Collect the t_ReturnType objects in t_ArgumentType
      //! @param PROTEIN_MODEL entity that contains a t_ReturnType
      //! @return returns Group of the collected t_ReturnType objects
      virtual util::SiPtrList< const ProteinModel> Collect( const ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class CollectorProteinModelConformations

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATIONS_H_ 
