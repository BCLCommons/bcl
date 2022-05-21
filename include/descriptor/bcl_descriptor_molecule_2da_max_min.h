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

#ifndef BCL_DESCRIPTOR_MOLECULE_2DA_MAX_MIN_H_
#define BCL_DESCRIPTOR_MOLECULE_2DA_MAX_MIN_H_

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
    //! @class Molecule2DAMaxMin
    //! @brief Provides methods for calculating 2D AutoCorrelation related input code for ANNs and SVMs
    //!
    //! @see @link example_descriptor_molecule_2da_max_min.cpp @endlink
    //! @author raftersa, mendenjl
    //! @date Jun 16, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule2DAMaxMin :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    public:

    /////////
    //enums//
    /////////

      enum Method
      {
        e_2DAMax, //! Whether to do max
        e_2DAMin  //! Whether to do min
      };

    private:

    //////////
    // data //
    //////////

      graph::ConstGraph< size_t, size_t>  m_Graph;        //!< Graph of the current molecule
      size_t                              m_NumberSteps;  //!< number of steps in 2DA function
      CheminfoProperty                    m_AtomProperty; //!< the atom property encode in the 2DA function
      bool                                m_MaxOrMin;     //!< tells functionality max or min
      float                               m_SubValue;     //!< provides the value to replace "empty" bins or
                                                          //!< tells the descriptor to use the worst seen value

    public:

      //! two instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_MaxInstance;
      static const util::SiPtr< const util::ObjectInterface> s_MinInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor, max vs min choice only
      Molecule2DAMaxMin( const Method MAXMIN);

      //! @brief constructor from number of steps, mapped atom property, and min/max choice
      Molecule2DAMaxMin
      (
        const size_t NUMBER_STEPS,
        const CheminfoProperty &ATOM_PROPERTY,
        const Method MAXMIN,
        const float SUBVALUE
      );

      //! @brief virtual copy constructor
      Molecule2DAMaxMin *Clone() const;

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
        return m_NumberSteps;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

      //! @brief get sub value of code
      //! @return sub value given to 2da code
      const float &GetSubValue() const;

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

    }; // class Molecule2DAMaxMin

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_2DA_MAX_MIN_H_
