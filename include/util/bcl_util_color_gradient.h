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

#ifndef BCL_UTIL_COLOR_GRADIENT_H_
#define BCL_UTIL_COLOR_GRADIENT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_function_interface_serializable.h"
#include "bcl_util_object_instances.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ColorGradient
    //! @brief ColorGradient class is for creating a taking a value and providing a corresponding color.
    //! @details It takes a double value (e.g. score) and returns the color point at which the value falls within the defined
    //! gradient based on min and max values. You provide an arbitrary number of color points which you want to
    //! transition between.
    //! The gradient between colors is linear. The range between the min and max values is split up evenly so that each
    //! color transition takes place over 1/(number color transitions) of the min to max value range.
    //! The color points are typically thought of as Red, Green, Blue format, although there
    //! is nothing in this class hard coded to force this. The color transitions take place in the order in which
    //! they are provided.
    //!
    //! The double is the value which will be converted into color points
    //! The color is an instance Colors enumerator
    //!
    //! @see @link example_util_color_gradient.cpp @endlink
    //! @author alexanns, karakam
    //! @date September 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ColorGradient :
      public FunctionInterfaceSerializable< double, linal::Vector3D>
    {
    protected:

    //////////
    // data //
    //////////

      //! range of the value that will be converted into a color
      math::Range< double> m_Range;

      //! All of the colors which will be used in the color gradient (N colors = N - 1 transitions)
      storage::Vector< Color> m_GradientColors;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ColorGradient();

      //! @brief constructor from a range and vector of gradient colors
      //! @param RANGE range of values
      //! @param GRADIENT_POINTS all of the colors which will be used in the color gradient
      ColorGradient
      (
        const math::Range< double> &RANGE,
        const storage::Vector< Color> &GRADIENT_POINTS
      );

      //! @brief virtual copy constructor
      //! @return new copy of this class
      ColorGradient *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a value and returning the corresponding color point
      //! @param VALUE the value which will be converted into a color
      //! @return returns the corresponding color
      linal::Vector3D operator()( const double &VALUE) const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // ColorGradient

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_COLOR_GRADIENT_H_
