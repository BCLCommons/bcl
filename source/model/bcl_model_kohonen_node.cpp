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
#include "model/bcl_model_kohonen_node.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> KohonenNode::s_Instance
    (
      GetObjectInstances().AddInstance( new KohonenNode())
    );

    //! @brief constructor from parameters
    //! @param POSITION the location of this node within the network.
    //! @param FEATURE_VECTOR the initial feature vector for this node
    //! @param RESULT_VECTOR the initial result vector for this node
    KohonenNode::KohonenNode
    (
      const linal::Vector< float> &POSITION,
      const linal::Vector< float> &FEATURE_VECTOR,
      const linal::Vector< float> &RESULT_VECTOR
    ) :
      m_Position( POSITION),
      m_FeatureAverage( FEATURE_VECTOR),
      m_ResultAverage( RESULT_VECTOR)
    {
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KohonenNode::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Position, ISTREAM);
      io::Serialize::Read( m_FeatureAverage, ISTREAM);
      io::Serialize::Read( m_ResultAverage, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KohonenNode::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Position, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FeatureAverage, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ResultAverage, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief adds the data from another node to the current node.
    //! @param OTHER the node to add from
    //! @return a reference to this node
    KohonenNode &KohonenNode::operator +=( const KohonenNode &OTHER)
    {
      if( m_FeatureAverage.GetWeight() == 0.0)
      {
        if( OTHER.m_FeatureAverage.GetWeight() > 0.0)
        {
          m_FeatureAverage = OTHER.m_FeatureAverage;
          m_ResultAverage = OTHER.m_ResultAverage;
        }
      }
      else if( OTHER.GetWeight() > 0.0)
      {
        m_FeatureAverage.AddWeightedObservation
        (
          OTHER.m_FeatureAverage.GetAverage(),
          OTHER.m_FeatureAverage.GetWeight()
        );
        m_ResultAverage.AddWeightedObservation
        (
          OTHER.m_ResultAverage.GetAverage(),
          OTHER.m_ResultAverage.GetWeight()
        );
      }
      return *this;
    }

    //! @brief adds the data from another node to the current node with a specified weight
    //! @param OTHER the node to add from
    //! @param WEIGHT the weight to given the reference vector from the newly added node
    //! @return a reference to this node
    KohonenNode &KohonenNode::AddNodeWithWeight( const KohonenNode &OTHER, const float &WEIGHT)
    {
      if( WEIGHT != 0.0)
      {
        m_FeatureAverage.AddWeightedObservation( OTHER.m_FeatureAverage.GetAverage(), WEIGHT);
        m_ResultAverage.AddWeightedObservation( OTHER.m_ResultAverage.GetAverage(), WEIGHT);
      }
      return *this;
    }

  } // namespace model
} // namespace bcl

