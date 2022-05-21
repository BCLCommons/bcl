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

#ifndef BCL_MATH_OBJECT_STOCHASTIC_SELECTOR_H_
#define BCL_MATH_OBJECT_STOCHASTIC_SELECTOR_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectStochasticSelector
    //! @brief ObjectStochasticSelector is a class that chooses between a set of objects, and is designed to be more
    //! @brief generalizable than the MutateDecisionNode class.  Note that any object that is used with this class
    //! @brief must be serializable.
    //!
    //! TODO: Combine this with ObjectProbabilityDistribution which is basically the same thing but not serializable
    //! @see @link example_math_object_stochastic_selector.cpp @endlink
    //! @author geanesar
    //! @date Sept 26 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ObjectStochasticSelector :
      public util::SerializableInterface
    {
    public:
    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

    //////////
    // data //
    //////////

      //! Map associating Implementations with their probabilities
      storage::Map< util::Implementation< t_DataType>, double> m_Distribution;

      //! Vector containing object data labels
      storage::Vector< util::Implementation< t_DataType> > m_DataLabels;

      //! @brief Vector containing probabilities
      storage::Vector< double> m_Probabilities;

      double m_Sum;

    private:

      //! @brief normalizes the probabilities of the map to 1.0
      void NormalizeProbabilities()
      {
        double sum( 0.0);

        for
        (
          typename storage::Map< util::Implementation< t_DataType>, double>::const_iterator itr_map( m_Distribution.Begin()),
            itr_map_end( m_Distribution.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          sum += itr_map->second;
        }

        double ratio( 1.0 / sum);

        for
        (
          typename storage::Map< util::Implementation< t_DataType>, double>::iterator itr_map( m_Distribution.Begin()),
            itr_map_end( m_Distribution.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          itr_map->second *= ratio;
        }
      }

      //! @brief multiplies all probabilities by a multiplier
      //! @param MULTIPLIER the multiplier to use
      void MultiplyProbabilities( const double &MULTIPLIER)
      {
        for
        (
          typename storage::Map< util::Implementation< t_DataType>, double>::iterator itr_map( m_Distribution.Begin()),
            itr_map_end( m_Distribution.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          itr_map->second *= MULTIPLIER;
        }
      }

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectStochasticSelector() :
        m_Distribution(),
        m_Sum( 0.0)
      {
      }

      //! @brief construct from a map 
      //! @param DISTRIBUTION ObjectProbabilityDistribution
      explicit ObjectStochasticSelector
      (
        const storage::Map< util::Implementation< t_DataType>, double> &DISTRIBUTION
      ) :
        m_Distribution( DISTRIBUTION),
        m_Sum( 0.0)
      {
      }

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      ObjectStochasticSelector< t_DataType> *Clone() const
      {
        return new ObjectStochasticSelector< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns class name when used in a dynamic context
      //! @return the class name when used in a dynamic context
      const std::string &GetAlias() const
      {
        static std::string s_alias( "StochasticSelector");
        return s_alias;
      }

      //! @brief gets the distribution used to select objects
      //! @return the map of objects (keys) to their probabilities (values)
      const storage::Map< util::Implementation< t_DataType>, double> &GetMap() const
      {
        return m_Distribution;
      }

      //! @brief get the number of available options
      //! @return the size of the decision node
      size_t GetSize() const
      {
        return m_Distribution.GetSize();
      }

    ////////////////
    // operations //
    ////////////////

    public:

      //! @brief add an object with probability
      //! @param OBJECT_LABEL the label of the implementation
      //! @param PROBABILITY the probability for that implementation to be selected
      void AddImplementation
      (
        const util::ObjectDataLabel &OBJECT_LABEL,
        const double &PROBABILITY
      )
      {
        util::Implementation< t_DataType> impl;
        impl->TryRead( OBJECT_LABEL);
        BCL_Assert( impl.IsDefined(), "Could not read label: " + util::Format()( OBJECT_LABEL));
        AddImplementation( impl, PROBABILITY);
      }

      //! @brief add an object with probability
      //! @param OBJECT_IMPL the implementation to add
      //! @param PROBABILITY the probability for that implementation to be selected
      void AddImplementation
      (
        const util::Implementation< t_DataType> &OBJECT_IMPL,
        const double &PROBABILITY
      )
      {
        MultiplyProbabilities( m_Sum / ( m_Sum + PROBABILITY));
        m_Distribution[ OBJECT_IMPL] = PROBABILITY / ( m_Sum + PROBABILITY);
        m_Sum += PROBABILITY;
      }

      //! @brief add an object with probability
      //! @param OBJECT the object to add to the map (cloned to an Implementation)
      //! @param PROBABILITY the probability for that implementation
      void AddImplementation
      (
        const t_DataType &OBJECT,
        const double &PROBABILITY
      )
      {
        AddImplementation( util::Implementation< t_DataType>( OBJECT), PROBABILITY);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief stochastically selects a member from the map based on relative probabilities
      //! @return the class that was selected
      const t_DataType &SelectRandomCase() const
      {
        BCL_Assert( m_Distribution.GetSize(), "There is nothing to select from.  Selector must have at least one member.");

        // random between 0-1.0
        double rand( random::GetGlobalRandom().Double());

        typename storage::Map< util::Implementation< t_DataType>, double>::const_iterator itr_map( m_Distribution.Begin()),
            itr_map_end( m_Distribution.End());

        // Equivalent to adding up probabilities until over the random number
        for( ; itr_map != itr_map_end; ++itr_map)
        {
          rand -= itr_map->second;
          if( rand <= 0)
          {
            break;
          }
        }
        return *( itr_map->first);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        if( m_DataLabels.GetSize() != m_Probabilities.GetSize())
        {
          BCL_MessageCrt
          (
            "ERROR: " + util::Format()( m_DataLabels.GetSize()) + " options and "
            + util::Format()( m_Probabilities.GetSize()) + " probabilities were specified.  Numbers do not match"
          );
          return false;
        }

        for( size_t num( 0); num < m_DataLabels.GetSize(); ++num)
        {
          AddImplementation( m_DataLabels( num), m_Probabilities( num));
        }

        // These vectors are no longer needed
        m_DataLabels.Resize( 0);
        m_Probabilities.Resize( 0);

        // Probabilities came in all at once, normalize them
        NormalizeProbabilities();
        return true;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {

        io::Serializer member_info;

        member_info.AddInitializer
        (
          "options",
          "mapping of objects to probabilities",
          io::Serialization::GetAgent( &m_DataLabels)
        );

        member_info.AddInitializer
        (
          "probabilities",
          "probability of each object in the \"objects\" list",
          io::Serialization::GetAgent( &m_Probabilities)
        );

        return member_info;
      }

    }; // template class ObjectStochasticSelector

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_OBJECT_STOCHASTIC_SELECTOR_H_
