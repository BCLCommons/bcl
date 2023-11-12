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

#ifndef BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_
#define BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreFunctionGeneric
    //! @brief Initializes a score function as a CheminfoProperty, computes that property on the given molecule, and
    //! returns one of several metrics from the property vector.
    //!
    //! @see @link example_chemistry_score_function_generic.cpp @endlink
    //! @author brownbp1
    //! @date May 01, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreFunctionGeneric :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {

    public:

    //////////
    // Enum //
    //////////

      //! @brief methods for how molecules should be selected
      enum CalculationType
      {
        e_Index = 0,               //!< Return a property value by a specific vector index
        e_Sum,                     //!< Return the sum of all property values
        e_Mean,                    //!< Return the average of all property values
        e_Min,                     //!< Return the minimum of all property values
        e_Max,                     //!< Return the maximum of all property values
        e_NormMax,                 //!< Return the normalized maximum of all property values
        e_SoftMax,                 //!< Return the exponentiated normalized maximum of all property values
        e_Entropy,                 //!< Return the information entropy computed from all values
        s_NumberCalculationTypes
      };

      //! @brief CalculationType as string
      //! @param CALCULATION_TYPE the calculation type whose name is desired
      //! @return the name as string
      static const std::string &GetCalculationTypeName( const CalculationType &CALCULATION_TYPE);

      //! ConjugationEnum simplifies the usage of the Conjugation enum of this class
      typedef util::WrapperEnum< CalculationType, &GetCalculationTypeName, s_NumberCalculationTypes> CalculationTypeEnum;

    private:

      //! the descriptor to use
      descriptor::CheminfoProperty m_Descriptor;

      //! the method to use when returning the descriptor value
      CalculationTypeEnum m_CalculationType;

      //! the reference index in the property vector for certain scores
      size_t m_PropertyIndex;

      //! invert the property values; occurs prior to any normalization
      bool m_Invert;

      //! normalize the property values; occurs after any inversion
      bool m_Normalize;

      //! add some noise to bins so that ln(0) is not undefined
      double m_Noise;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScoreFunctionGeneric();

      //! @brief constructor with a property parameter
      //! @param DESCRIPTOR the descriptor to use
      ScoreFunctionGeneric
      (
        const descriptor::CheminfoProperty &DESCRIPTOR
      );

      //! @brief constructor with all parameters
      //! @param DESCRIPTOR the descriptor to use
      //! @param INDEX the reference index for certain scores
      //! @param INVERT invert each value in the property array
      //! @param NORMALIZE normalize property array to sum of values
      //! @param NOISE add some small value to bins to avoid ln(0) = nan
      explicit ScoreFunctionGeneric
      (
        const descriptor::CheminfoProperty &DESCRIPTOR,
        const size_t INDEX,
        const bool INVERT,
        const bool NORMALIZE,
        const double NOISE
      );

      //! @brief Clone function
      //! @return pointer to new ScoreFunctionGeneric
      ScoreFunctionGeneric *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the class name when used in a dynamic context
      //! @return the class name when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return the value at a single index
      const double CalcIndexValue( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the sum of property values
      const double CalcSum( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the mean descriptor value
      const double CalcMean( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the min descriptor value
      const double CalcMin( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the max descriptor value
      const double CalcMax( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the maximum value after prop[ref_index]/sum_0-->N(prop)
      const double CalcNormMax( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return exp(prop[ref_index])/sum_0-->N(exp(prop))
      const double CalcSoftMax( const linal::Vector< float> &PROPERTIES) const;

      //! @brief return the entropy of the dataset
      const double CalcEntropy( const linal::Vector< float> &PROPERTIES) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return propensity score for observing rotamers that exist in conformation
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief invert each value of the property vector
      void Invert( linal::Vector< float> PROPERTIES) const;

      //! @brief normalize property vector by sum of all values
      void Normalize( linal::Vector< float> PROPERTIES) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief gets a serializer for constructing this class in a dynamic context
      //! @return the serializer containing member data
      io::Serializer GetSerializer() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ScoreFunctionGeneric

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_
