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

#ifndef BCL_DESCRIPTOR_FOR_EACH_H_
#define BCL_DESCRIPTOR_FOR_EACH_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ForEach
    //! @brief Substitutes each of the given values into a label to create a new descriptor
    //!
    //! @see @link example_descriptor_for_each.cpp @endlink
    //! @author mendenjl
    //! @date Feb 02, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class ForEach :
      public Combine< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      //! Template label
      util::ObjectDataLabel m_TemplateLabel;

      //! Variable name for the template
      std::string m_Variable;

      //! Values that the variable can take on
      storage::Vector< util::ObjectDataLabel> m_Values;

      //! Value descriptor(s). Allows for use of previously defined lists, e.g. X in Define(X=Combine(a1,a2...))
      //! The implementation must be dynamically-castable to a Combine (list)
      storage::Vector< util::ObjectDataLabel> m_ValuesDescriptors;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      ForEach *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief interpret the variables in the given combine as arguments to be used in the template
      //! @param BASE the base to use. If it is a defined list or list directly, the internal descriptors will be used instead
      //! @param ERR_STREAM stream to write out errors to
      bool AddDescriptorArguments
      (
        const Base< t_DataType, t_ReturnType> &BASE,
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    }; // class ForEach

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API ForEach< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_FOR_EACH_H_
