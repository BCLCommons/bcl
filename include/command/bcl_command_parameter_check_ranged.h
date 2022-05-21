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

#ifndef BCL_COMMAND_PARAMETER_CHECK_RANGED_H_
#define BCL_COMMAND_PARAMETER_CHECK_RANGED_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckRanged
    //! @brief This class is a command line helper template class derived from ParameterCheck that checks, whether a
    //! given parameter is within a certain range. You may use every template that has a definition for < and ==.
    //!
    //! @see @link example_command_parameter_check_ranged.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ParameterCheckRanged :
      public ParameterCheckInterface
    {

    private:

    //////////
    // data //
    //////////

      math::Range< t_DataType> m_Range; //!< contains minimum and maximum allowed parameter

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ParameterCheckRanged()
      {
        const t_DataType min( std::numeric_limits< t_DataType>::min()), max( std::numeric_limits< t_DataType>::max());
        const bool is_signed( std::numeric_limits< t_DataType>::is_signed);

        // if t_Datatype is signed, the min is -max
        m_Range = is_signed ? math::Range< t_DataType>( -max, max) : math::Range< t_DataType>( min, max);
      }

      //! @brief construct from allowed parameters
      //! @param MIN minimum value to be allowed
      //! @param MAX maximum value to be allowed
      ParameterCheckRanged( const t_DataType MIN, const t_DataType MAX)
      {
        m_Range = math::Range< t_DataType>( MIN, MAX);
      }

      //! @brief construct from allowed parameters
      //! @param MIN minimum value to be allowed
      explicit ParameterCheckRanged( const t_DataType MIN)
      {
        m_Range = math::Range< t_DataType>( MIN, std::numeric_limits< t_DataType>::max());
      }

      //! @brief virtual copy constructor
      ParameterCheckRanged< t_DataType> *Clone() const
      {
        return new ParameterCheckRanged< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return const ref std::string - the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is allowed, i.e. in the specified range
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      bool IsAllowedParameter
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const
      {
        if( util::IsNumerical( PARAMETER))
        {
          const t_DataType numerical_value( util::ConvertStringToNumericalValue< t_DataType>( PARAMETER));

          // return true if numerical_value is in range m_Min..m_Max
          if( m_Range.IsWithin( numerical_value))
          {
            return true;
          }
          // else write to the error message
          else
          {
            ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is not in the range " << m_Range.GetString() << '\n';
          }
        }
        else
        {
          ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is not a numerical value" << '\n';
        }

        return false;
      } // IsAllowedParameter

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM - the outstream to write to
      //! @param INDENT amount to indent each line after the first
      //! @return return the stream after you write to it
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        // write range
        OSTREAM << "range: " << m_Range.GetString();

        // end
        return OSTREAM;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Range, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      //! @see ParameterCheckInterface::Write
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Range, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class ParameterRanged

    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> ParameterCheckRanged< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckRanged< t_DataType>())
    );

    BCL_EXPIMP_TEMPLATE template class BCL_API ParameterCheckRanged< double>;

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_RANGED_H_
