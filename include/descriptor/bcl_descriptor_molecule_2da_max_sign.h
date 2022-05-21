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

#ifndef BCL_DESCRIPTOR_MOLECULE_2DA_MAX_SIGN_H_
#define BCL_DESCRIPTOR_MOLECULE_2DA_MAX_SIGN_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule2DAMaxSign
    //! @brief Provides methods for calculating 2D AutoCorrelation related input code for ANNs and SVMs
    //!
    //! @see @link example_descriptor_molecule_2da_max_sign.cpp @endlink
    //! @author raftersa, mendenjl
    //! @date Jun 23, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule2DAMaxSign :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      graph::ConstGraph< size_t, size_t>  m_Graph;         //!< Graph of the current molecule
      CheminfoProperty                    m_AtomProperty;  //!< the atom property encode in the 3D autocorrelation function
      size_t                              m_NumberSteps;   //!< number of steps in 2D autocorrelation function

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule2DAMaxSign();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule2DAMaxSign
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t NUMBER_STEPS
      );

      //! @brief virtual copy constructor
      Molecule2DAMaxSign *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_NumberSteps * 3;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    protected:

      //! @brief add an observed distance/property value to m_DiscreteCode
      //! @param DISTANCE actual distance of the two atoms
      //! @param PROP_A property from atom A
      //! @param PROP_B property from atom B
      void Accumulate
      (
        const size_t &DISTANCE,
        const float &PROP_A,
        const float &PROP_B
      );

    }; // class Molecule2DAMaxSign

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_2DA_MAX_SIGN_H_
