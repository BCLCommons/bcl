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

// include header of this class
#include "util/bcl_util_color_gradient.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const SiPtr< const ObjectInterface> ColorGradient::s_Instance
    (
      Enumerated< FunctionInterfaceSerializable< double, linal::Vector3D> >::AddInstance( new ColorGradient())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ColorGradient::ColorGradient() :
      m_Range( 0.0, 1.0),
      m_GradientColors( GetColors().GetRainbow())
    {
    }

    //! @brief constructor from a range and vector of gradient colors
    //! @param RANGE range of values
    //! @param GRADIENT_POINTS all of the colors which will be used in the color gradient
    ColorGradient::ColorGradient
    (
      const math::Range< double> &RANGE,
      const storage::Vector< Color> &GRADIENT_POINTS
    ) :
      m_Range( RANGE),
      m_GradientColors( GRADIENT_POINTS)
    {
      BCL_Assert( GRADIENT_POINTS.GetSize() > 1, "The gradient needs to have more than one color");
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    ColorGradient *ColorGradient::Clone() const
    {
      return new ColorGradient( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ColorGradient::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ColorGradient::GetAlias() const
    {
      static const std::string s_Name( "ColorGradient");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking a value and returning the corresponding color point
    //! @param VALUE the value which will be converted into a color
    //! @return returns the corresponding color
    linal::Vector3D ColorGradient::operator()( const double &VALUE) const
    {
      // if value larger than the upper bound
      if( VALUE >= m_Range.GetMax())
      {
        BCL_MessageCrt
        (
          Format()( VALUE) + " is not within range " + m_Range.GetString() + " returning max color"
        );

        // return the color for the largest expected value
        return m_GradientColors.LastElement();
      }
      // if value is smaller than the lower bound
      else if( VALUE < m_Range.GetMin())
      {
        BCL_MessageCrt
        (
          Format()( VALUE) + " is not within range " + m_Range.GetString() + " returning min color"
        );

        // return the color for the smallest expected value
        return m_GradientColors.FirstElement();
      }

      // rescale the value according to a new range
      const double position( ( VALUE - m_Range.GetMin()) / m_Range.GetWidth() * ( m_GradientColors.GetSize() - 1));

      // calculate the index of the color in the gradient this VALUE corresponds to
      const size_t interval_index( ( size_t)( position));

      // calculate the difference between start and end points of this gradient
      const linal::Vector3D color_difference( *m_GradientColors( interval_index + 1) - *m_GradientColors( interval_index));

      // calculate the new color and return it
      return *m_GradientColors( interval_index) + color_difference * ( position - double( interval_index));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ColorGradient::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Range, ISTREAM);
      io::Serialize::Read( m_GradientColors, ISTREAM);
      BCL_Assert( m_GradientColors.GetSize() > 1, "The gradient needs at least 2 colors!");

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ColorGradient::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Range, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GradientColors, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ColorGradient::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "For creating a taking a value and providing a corresponding color."
        "It takes a double value (e.g. score) and returns the color point at which the value falls within the defined"
        "gradient based on min and max values. You provide an arbitrary number of color points which you want to"
        "transition between."
        "The gradient between colors is linear. The range between the min and max values is split up evenly so that each"
        "color transition takes place over 1/(number color transitions) of the min to max value range."
        "The color points are typically thought of as Red, Green, Blue format, although there"
        "is nothing in this class hard coded to force this. The color transitions take place in the order in which"
        "they are provided. The double is the value which will be converted into color points."
        "The color is an instance Colors enumerator."
      );

      member_data.AddInitializer
      (
        "range",
        "the range of values over which a color code will be applied (must be quoted if range string contains () or ,)",
        io::Serialization::GetAgent( &m_Range),
        "\"(0.0,0.0)\""
      );

      member_data.AddInitializer
      (
        "color_list",
        "The list of colors that the color gradient will use. First color is for range minima, last is for range max.",
        io::Serialization::GetAgentWithSizeLimits( &m_GradientColors, 2),
        ObjectDataLabel
        (
          storage::Vector< ObjectDataLabel>::Create
          (
            ObjectDataLabel( GetColors().e_Black->GetName()),
            ObjectDataLabel( GetColors().e_White->GetName())
          )
        ).ToString()
      );

      return member_data;
    }

  } // namespace util
} // namespace bcl

