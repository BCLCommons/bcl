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
#include "model/bcl_model_feature_label_set.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FeatureLabelSet::s_Instance
    (
      GetObjectInstances().AddInstance( new FeatureLabelSet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FeatureLabelSet::FeatureLabelSet() :
      m_OrderedProperties(),
      m_Size( 0),
      m_OuterName( "Combine")
    {
    }

    //! @brief constructor from outer name
    //! @param NAME name for the set, usually a descriptor that takes a list of other descriptors
    FeatureLabelSet::FeatureLabelSet
    (
      const std::string &NAME,
      const util::Implementation< util::ImplementationInterface> &IMPL
    ) :
      m_OrderedProperties(),
      m_Size( 0),
      m_OuterName( NAME),
      m_ImplementationInterface
      (
        IMPL.IsDefined()
        ? util::Implementation< util::ImplementationInterface>( IMPL->Empty())
        : IMPL
      )
    {
      if( NAME.find( '!') != std::string::npos && !m_ImplementationInterface.IsDefined())
      {
        storage::Vector< std::string> split( util::SplitString( NAME, "!"));
        if( split.GetSize() > 1 && !split( 0).empty())
        {
          std::stringstream err_stream;
          const std::string impl_interface( split( 0));
          split.RemoveElements( 0, 1);
          if( m_ImplementationInterface.TryRead( util::ObjectDataLabel( impl_interface), err_stream))
          {
            m_OuterName = util::Join( "!", split);
          }
        }
        else
        {
          m_OuterName = split( 0);
        }
      }
    }

    //! @brief constructor from members
    //! @param NAME name for the set, usually a descriptor that takes a list of other descriptors
    //! @param PROPERTIES properties of the feature
    //! @param PROPERTY_SIZES sizes of those properties
    FeatureLabelSet::FeatureLabelSet
    (
      const std::string &NAME,
      const storage::Vector< util::ObjectDataLabel> &PROPERTIES,
      const storage::Vector< size_t> &PROPERTY_SIZES,
      const util::Implementation< util::ImplementationInterface> &IMPL
    ) :
      m_OrderedProperties(),
      m_Size( 0),
      m_OuterName( NAME),
      m_ImplementationInterface
      (
        IMPL.IsDefined()
        ? util::Implementation< util::ImplementationInterface>( IMPL->Empty())
        : IMPL
      )
    {
      BCL_Assert
      (
        PROPERTIES.GetSize() == PROPERTY_SIZES.GetSize(),
        "Must have a size for each property, but received "
        + util::Format()( PROPERTIES.GetSize()) + " properties and " + util::Format()( PROPERTY_SIZES.GetSize())
        + " sizes: " + util::Format()( PROPERTIES) + " " + util::Format()( PROPERTY_SIZES)
      );
      for
      (
        size_t property_id( 0), number_properties( PROPERTIES.GetSize());
        property_id < number_properties;
        ++property_id
      )
      {
        PushBack( PROPERTIES( property_id), PROPERTY_SIZES( property_id));
      }
      if( NAME.find( '!') != std::string::npos && !m_ImplementationInterface.IsDefined())
      {
        storage::Vector< std::string> split( util::SplitString( NAME, "!"));
        if( split.GetSize() > 1 && !split( 0).empty())
        {
          std::stringstream err_stream;
          const std::string impl_interface( split( 0));
          split.RemoveElements( 0, 1);
          if( m_ImplementationInterface.TryRead( util::ObjectDataLabel( impl_interface), err_stream))
          {
            m_OuterName = util::Join( "!", split);
          }
        }
        else
        {
          m_OuterName = split( 0);
        }
      }
    }

    //! @brief virtual copy constructor
    FeatureLabelSet *FeatureLabelSet::Clone() const
    {
      return new FeatureLabelSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FeatureLabelSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the alias / name used when constructing this object
    //! @returnthe alias / name used when constructing this object
    const std::string &FeatureLabelSet::GetAlias() const
    {
      return m_OuterName;
    }

    //! @brief return the data label
    //! @return data label as string
    std::string FeatureLabelSet::GetString() const
    {
      return GetLabel().ToString();
    }

    //! @brief return a label with the descriptors in this
    //! @return a label with the descriptors in this
    util::ObjectDataLabel FeatureLabelSet::GetLabel() const
    {
      if( m_OuterName.empty() && m_OrderedProperties.GetSize() == size_t( 1))
      {
        // return the only descriptors label
        return m_OrderedProperties.FirstElement();
      }
      return util::ObjectDataLabel( m_OuterName, m_OrderedProperties);
    }

    //! @brief Get a vector containing all member data from this property as separate data labels
    const storage::Vector< util::ObjectDataLabel> &FeatureLabelSet::GetMemberLabels() const
    {
      return m_OrderedProperties;
    }

    //! @brief get the sizes of each property
    //! @return sizes of each property, in the same order as GetMemberLabels
    const storage::Vector< size_t> &FeatureLabelSet::GetPropertySizes() const
    {
      return m_PropertiesSizes;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief push back a property and its expected size
    //! @param PROPERTY_LABEL label describing the property
    //! @param PROPERTY_SIZE size of the label
    void FeatureLabelSet::PushBack
    (
      const util::ObjectDataLabel &PROPERTY_LABEL,
      const size_t &PROPERTY_SIZE
    )
    {
      // remove the split label set, if it is present
      if( m_SplitFeatures.IsDefined())
      {
        m_SplitFeatures = util::OwnPtr< FeatureLabelSet>();
      }

      // get the range for the last element in ordered properties, if applicable

      // add the property to the vector
      m_OrderedProperties.PushBack( PROPERTY_LABEL);

      // create a math::Range object to hold the range of this property
      math::Range< size_t> range
      (
        math::RangeBorders::e_LeftClosed, m_Size, m_Size + PROPERTY_SIZE, math::RangeBorders::e_RightClosed
      );

      // add the size to the total size
      m_Size += PROPERTY_SIZE;

      // add the size to the property size vector
      m_PropertiesSizes.PushBack( PROPERTY_SIZE);

      // add the property to the map
      m_PropertiesToRanges[ PROPERTY_LABEL] = range;

      BCL_Assert
      (
        m_OrderedProperties.GetSize() == m_PropertiesToRanges.GetSize(),
        "Duplicated object data labels are not allowed in feature label set. Duplicate label was: "
        + PROPERTY_LABEL.ToString()
      );
    }

    //! @param LABEL label of the desired property
    //! @return segment(s) of feature vector used by the given property
    storage::Vector< size_t> FeatureLabelSet::GetPropertyIndices( const util::ObjectDataLabel &LABEL) const
    {
      // look for the label in the map
      storage::Map< util::ObjectDataLabel, math::Range< size_t> >::const_iterator itr( m_PropertiesToRanges.Find( LABEL));

      // if the property is not found, consider whether it is a subproperty
      if( itr == m_PropertiesToRanges.End())
      {
        // property did not exist directly; try to retrieve it assuming it is a subproperty
        return LABEL == GetLabel() ? storage::CreateIndexVector( m_Size) : GetSubPropertyIndices( LABEL);
      }

      // return indices in the selected range
      return storage::CreateIndexVector( itr->second.GetMin(), itr->second.GetMax());
    }

    //! @brief create a feature label subset
    //! @param SUB_FEATURES indices of the complete feature to keep
    //! @return feature label set with the given sub features
    //! @note: if a descriptor has multiple return values and must be split, it is assumed that there is a
    //! meta-descriptor named Partial, that takes another descriptor and a parameter called indices, which contains
    //! the indices of the feature labels to keep.
    FeatureLabelSet FeatureLabelSet::CreateSubFeatureLabelSet( const storage::Vector< size_t> &SUB_FEATURES) const
    {
      FeatureLabelSet desired_features( m_OuterName, m_ImplementationInterface);

      // determine which features are needed
      storage::Vector< size_t> index_in_kept_features( m_Size, util::GetUndefined< size_t>());
      size_t counter( 0);
      for
      (
        storage::Vector< size_t>::const_iterator
          itr( SUB_FEATURES.Begin()), itr_end( SUB_FEATURES.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        index_in_kept_features( *itr) = counter;
      }

      storage::Vector< util::ObjectDataLabel> new_labels( m_PropertiesToRanges.GetSize());
      storage::Vector< size_t> new_sizes( m_PropertiesToRanges.GetSize(), size_t( 0));

      storage::Map< util::ObjectDataLabel, size_t> label_to_original_index;
      counter = 0;
      for
      (
        util::ObjectDataLabel::const_iterator
          itr( m_OrderedProperties.Begin()), itr_end( m_OrderedProperties.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        label_to_original_index[ *itr] = counter;
        // always add zero-sized properties
        if( m_PropertiesSizes( counter) == size_t( 0))
        {
          desired_features.PushBack( *itr, size_t( 0));
        }
      }

      // iterate over all feature
      for
      (
        storage::Map< util::ObjectDataLabel, math::Range< size_t> >::const_iterator
          itr( m_PropertiesToRanges.Begin()), itr_end( m_PropertiesToRanges.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the size of this descriptor
        const size_t feature_size( itr->second.GetWidth());
        const size_t feature_start( itr->second.GetMin());

        // get the index of this feature in the original feature label set
        const size_t original_index( label_to_original_index[ itr->first]);

        storage::List< size_t> sub_features;
        for( size_t index( feature_start), index_end( itr->second.GetMax()); index < index_end; ++index)
        {
          if( util::IsDefined( index_in_kept_features( index)))
          {
            sub_features.PushBack( index - feature_start);
          }
        }

        // if part of this descriptor is to be kept, add it to the lists
        if( sub_features.GetSize())
        {
          new_sizes( original_index) = sub_features.GetSize();
          if( sub_features.GetSize() == feature_size)
          {
            // add the entire property to the desired features set
            new_labels( original_index) = CollapsePartials( itr->first);
          }
          else
          {
            // create the sub-feature label
            util::ObjectDataLabel subfeature
            (
              "Partial",
              storage::Vector< util::ObjectDataLabel>::Create
              (
                itr->first,
                util::ObjectDataLabel
                (
                  "indices",
                  io::Serialization::GetAgent( &sub_features)->GetLabel()
                )
              )
            );
            new_labels( original_index) = CollapsePartials( subfeature);
          }
        }
      }

      for( size_t label_index( 0), labels_size( new_sizes.GetSize()); label_index < labels_size; ++label_index)
      {
        if( new_sizes( label_index))
        {
          desired_features.PushBack( new_labels( label_index), new_sizes( label_index));
        }
      }

      return desired_features;
    }

    //! @brief create a feature label subset
    //! @param SKIP_ZERO_LENGTH_DESCRIPTORS true to ignore zero-length descriptors during the split. In this case, all
    //!        descriptors in the returned labels will have size == 1
    //! @return feature label set with each multi-column descriptor split using the Partial descriptor (vector index)
    //! @note: if a descriptor has multiple return values and must be split, it is assumed that there is a
    //! meta-descriptor named Partial, that takes another descriptor and a parameter called indices, which contains
    //! the indices of the feature labels to keep.
    FeatureLabelSet FeatureLabelSet::SplitFeatureLabelSet( const bool &SKIP_ZERO_LENGTH_DESCRIPTORS) const
    {
      const size_t number_zero_length_descriptors
      (
        SKIP_ZERO_LENGTH_DESCRIPTORS
        ? 0
        : std::count( m_PropertiesSizes.Begin(), m_PropertiesSizes.End(), size_t( 0))
      );
      storage::Vector< util::ObjectDataLabel> labels;
      labels.AllocateMemory( m_Size + number_zero_length_descriptors);

      storage::Vector< size_t> label_sizes( m_Size + number_zero_length_descriptors, size_t( 1));

      storage::Vector< size_t> sub_feature( 1, 0);
      util::OwnPtr< io::SerializationInterface> indices_container_label_sp
      (
        io::Serialization::GetAgent( &sub_feature)
      );

      storage::Vector< size_t>::const_iterator itr_sizes( m_PropertiesSizes.Begin());
      storage::Vector< size_t>::iterator itr_split_sizes( label_sizes.Begin());
      for
      (
        util::ObjectDataLabel::const_iterator
          itr( m_OrderedProperties.Begin()), itr_end( m_OrderedProperties.End());
        itr != itr_end;
        ++itr, ++itr_sizes, ++itr_split_sizes
      )
      {
        const size_t size( *itr_sizes);

        // handle single-value properties trivially
        if( size <= size_t( 1))
        {
          if( size == size_t( 0))
          {
            if( SKIP_ZERO_LENGTH_DESCRIPTORS)
            {
              --itr_split_sizes;
              continue;
            }
            else
            {
              *itr_split_sizes = 0;
            }
          }
          labels.PushBack( CollapsePartials( *itr));
          continue;
        }

        for( size_t id( 0); id < size; ++id)
        {
          sub_feature( 0) = id;
          // create the sub-feature label
          util::ObjectDataLabel subfeature
          (
            "Partial",
            storage::Vector< util::ObjectDataLabel>::Create
            (
              *itr,
              util::ObjectDataLabel( "indices", indices_container_label_sp->GetLabel())
            )
          );

          labels.PushBack( CollapsePartials( subfeature));
        }
      }

      return FeatureLabelSet( m_OuterName, labels, label_sizes, m_ImplementationInterface);
    }

    //! @brief get common overlap of features in this and given FeatureLabelSets
    FeatureLabelSet FeatureLabelSet::GetCommonFeatures( const FeatureLabelSet &COMPARE) const
    {
      std::set< util::ObjectDataLabel> labels_this
      (
        m_OrderedProperties.Begin(),
        m_OrderedProperties.End()
      );

      std::set< util::ObjectDataLabel> labels_that
      (
        COMPARE.m_OrderedProperties.Begin(),
        COMPARE.m_OrderedProperties.End()
      );

      std::list< util::ObjectDataLabel> common_labels_list;

      std::set_intersection
      (
        labels_this.begin(),
        labels_this.end(),
        labels_that.begin(),
        labels_that.end(),
        std::inserter( common_labels_list, common_labels_list.begin())
      );

      storage::Vector< util::ObjectDataLabel> common_labels( common_labels_list.begin(), common_labels_list.end());
      storage::Vector< size_t> common_sizes;
      for
      (
        util::ObjectDataLabel::const_iterator
          itr( common_labels.Begin()), itr_end( common_labels.End());
        itr != itr_end;
        ++itr
      )
      {
        common_sizes.PushBack( m_PropertiesToRanges.Find( *itr)->second.GetWidth());
      }
      return FeatureLabelSet( m_OuterName, common_labels, common_sizes, m_ImplementationInterface);
    }

    //! @brief merge two object data labels, respecting split descriptors
    //! @param LABEL_A, LABEL_B the object data labels of interest
    //! @return the merged labels
    util::ObjectDataLabel FeatureLabelSet::MergeConsideringPartials
    (
      const util::ObjectDataLabel &LABEL_A,
      const util::ObjectDataLabel &LABEL_B
    )
    {
      // for each descriptor that is a partial (e.g. column of another descriptor):
      // add the entry into the map, along with the corresponding descriptor indices
      // for each descriptor that is not a partial,
      // add the entry into the map without any indices
      storage::Map< util::ObjectDataLabel, storage::Set< size_t> > descriptor_indices;
      storage::Set< util::ObjectDataLabel> complete_descriptors;
      for
      (
        storage::Vector< util::ObjectDataLabel>::const_iterator itr_a( LABEL_A.Begin()), itr_a_end( LABEL_A.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        // attempt to decompose the descriptor into a partial
        storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> > decomposed( DecomposePartial( *itr_a));
        if( decomposed.First().IsEmpty())
        {
          // complete descriptor, add it to the set
          complete_descriptors.Insert( *itr_a);
          // if the descriptor was in the map, remove it
          descriptor_indices.Erase( *itr_a);
        }
        else if( !complete_descriptors.Contains( decomposed.First()))
        {
          const storage::Vector< size_t> &descriptor_ids( decomposed.Second());
          // partial descriptor, add the indices to the map
          descriptor_indices[ decomposed.First()].InsertElements( descriptor_ids.Begin(), descriptor_ids.End());
        }
      }

      for
      (
        storage::Vector< util::ObjectDataLabel>::const_iterator itr_b( LABEL_B.Begin()), itr_b_end( LABEL_B.End());
        itr_b != itr_b_end;
        ++itr_b
      )
      {
        // attempt to decompose the descriptor into a partial
        storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> > decomposed( DecomposePartial( *itr_b));
        if( decomposed.First().IsEmpty())
        {
          // complete descriptor, add it to the set
          complete_descriptors.Insert( *itr_b);
          // if the descriptor was in the map, remove it
          descriptor_indices.Erase( *itr_b);
        }
        else if( !complete_descriptors.Contains( decomposed.First()))
        {
          const storage::Vector< size_t> &descriptor_ids( decomposed.Second());
          // partial descriptor, add the indices to the map
          descriptor_indices[ decomposed.First()].InsertElements( descriptor_ids.Begin(), descriptor_ids.End());
        }
      }

      // get all the complete descriptors
      storage::Vector< util::ObjectDataLabel> descriptors( complete_descriptors.Begin(), complete_descriptors.End());

      storage::Vector< size_t> sub_feature( 1, 0);
      util::OwnPtr< io::SerializationInterface> indices_container_label_sp
      (
        io::Serialization::GetAgent( &sub_feature)
      );

      // add all other descriptors from the map
      for
      (
        storage::Map< util::ObjectDataLabel, storage::Set< size_t> >::const_iterator
          itr( descriptor_indices.Begin()), itr_end( descriptor_indices.End());
        itr != itr_end;
        ++itr
      )
      {
        // copy the indices
        sub_feature = storage::Vector< size_t>( itr->second.Begin(), itr->second.End());

        // compose the partial descriptor
        util::ObjectDataLabel subfeature
        (
          "Partial",
          storage::Vector< util::ObjectDataLabel>::Create
          (
            itr->first,
            util::ObjectDataLabel( "indices", indices_container_label_sp->GetLabel())
          )
        );
        descriptors.PushBack( CollapsePartials( subfeature));
      }
      return util::ObjectDataLabel( LABEL_A.GetName(), LABEL_A.GetValue(), descriptors);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FeatureLabelSet::Read( std::istream &ISTREAM)
    {
      // load the label
      util::ObjectDataLabel label( ISTREAM);

      // load the sizes of each property
      storage::Vector< size_t> property_sizes;
      ISTREAM >> property_sizes;

      m_ImplementationInterface.Reset();
      if( label.GetValue().find( '!') == std::string::npos)
      {
        // recreate this object with the properties and their sizes
        *this = FeatureLabelSet( label.GetValue(), label.GetArguments(), property_sizes);
      }
      else
      {
        storage::Vector< std::string> split( util::SplitString( label.GetValue(), "!"));
        std::stringstream err_stream;
        const std::string impl_interface( split( 0));
        split.RemoveElements( 0, 1);
        if( m_ImplementationInterface.TryRead( util::ObjectDataLabel( impl_interface), err_stream))
        {
          // valid object type
          *this =
            FeatureLabelSet
            (
              util::Join( "!", split),
              label.GetArguments(),
              property_sizes,
              m_ImplementationInterface
            );
        }
        else
        {
          // invalid object type
          BCL_MessageStd
          (
            "Warning; descriptor namespace with label: " + impl_interface
            + " could not be read; selection of descriptor subsets will not work if any descriptors had implicit / "
            "default values given. error output: " + err_stream.str()
          );
          *this = FeatureLabelSet( label.GetValue(), label.GetArguments(), property_sizes);
        }
      }

      // remove the split label set, if it is present
      if( m_SplitFeatures.IsDefined())
      {
        m_SplitFeatures = util::OwnPtr< FeatureLabelSet>();
      }

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FeatureLabelSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      OSTREAM <<
        util::ObjectDataLabel
        (
          ( m_ImplementationInterface.IsDefined() ? m_ImplementationInterface->GetAlias() : "") + "!" + m_OuterName,
          m_OrderedProperties
        ).ToString( 120, INDENT, 1) << '\n';
      OSTREAM << m_PropertiesSizes;

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief helper function to decompose a partial into the innermost object data label and the vector of indices
    //! @param LABEL the partial label of interest
    //! @return a pair containing the main descriptors object data label and the indices for the partial
    storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >
      FeatureLabelSet::DecomposePartial( const util::ObjectDataLabel &LABEL)
    {
      static const std::string s_partial( "Partial"), s_indices( "indices");

      // indices property, once found
      util::ObjectDataLabel::const_iterator itr_indices_property( LABEL.End());

      // storage for ptr to main property, once found
      util::ObjectDataLabel::const_iterator itr_main_property( LABEL.End());

      // first, test whether the property is a subproperty; which must have 2 arguments
      if( LABEL.GetValue() == s_partial && LABEL.GetNumberArguments() == size_t( 2))
      {
        itr_indices_property = LABEL.FindValue( s_indices, false);
        itr_main_property = LABEL.Begin();
        if( itr_indices_property == LABEL.Begin())
        {
          ++itr_main_property;
        }
      }

      // label was not a partial
      if( itr_indices_property == LABEL.End())
      {
        return storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >();
      }

      storage::Vector< size_t> indices_main_property;
      util::OwnPtr< io::SerializationInterface> indices_container_label_sp
      (
        io::Serialization::GetAgent( &indices_main_property)
      );

      std::stringstream err_stream;
      BCL_Assert
      (
        indices_container_label_sp->TryRead( *itr_indices_property, err_stream),
        "Non-numeric value in sub-property indices: " + err_stream.str()
      );
      return
        storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >( *itr_main_property, indices_main_property);
    }

    namespace
    {
      //! @brief helper function to determine synonyms for a given object data label by trying to initialize a descriptor
      //!        from the given label
      //! @param LABEL the label of interest
      //! @param SYNONYMS_CONTAINER container in which to store synonyms that are different from LABEL
      template< typename t_ObjectType, typename t_ReturnType>
      void TryInterpretingAsDescriptor
      (
        const util::ObjectDataLabel &LABEL,
        storage::Set< util::ObjectDataLabel> &SYNONYMS_CONTAINER
      )
      {
        std::ostringstream err_stream;
        util::Implementation< descriptor::Base< t_ObjectType, t_ReturnType> > impl( LABEL, err_stream);
        if( impl.IsDefined())
        {
          SYNONYMS_CONTAINER.Insert( LABEL);
        }
      }
    }

    //! @brief get a vector with all known synonyms for a given label, created from trying to update the data label
    //!        with all template instances of descriptor::Base
    //! @param LABEL the label to retrieve all synonyms for. These are primarily necessary when
    storage::Set< util::ObjectDataLabel> FeatureLabelSet::GetSynonyms( const util::ObjectDataLabel &LABEL)
    {
      storage::Set< util::ObjectDataLabel> synonyms;
      TryInterpretingAsDescriptor< chemistry::AtomConformationalInterface, char>( LABEL, synonyms);
      TryInterpretingAsDescriptor< chemistry::AtomConformationalInterface, float>( LABEL, synonyms);
      TryInterpretingAsDescriptor< biol::AABase, char>( LABEL, synonyms);
      TryInterpretingAsDescriptor< biol::AABase, float>( LABEL, synonyms);
      TryInterpretingAsDescriptor< char, char>( LABEL, synonyms);
      TryInterpretingAsDescriptor< char, float>( LABEL, synonyms);
      synonyms.Erase( LABEL);
      return synonyms;
    }

    //! @brief clean a label that might have multiple levels of partials
    //! @param LABEL the old label to consider
    util::ObjectDataLabel FeatureLabelSet::CollapsePartials( const util::ObjectDataLabel &LABEL)
    {
      static const std::string s_partial( "Partial"), s_indices( "indices");
      if( LABEL.GetValue() != s_partial || LABEL.GetNumberArguments() != size_t( 2))
      {
        return LABEL;
      }

      // decompose the partial into partial and indices
      storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >
        decomposed_partial_outer( DecomposePartial( LABEL));

      // determine whether there was a partial
      if( decomposed_partial_outer.First().IsEmpty() || decomposed_partial_outer.Second().IsEmpty())
      {
        return LABEL;
      }

      // find the sub-partials indices
      storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >
        decomposed_partial_inner( DecomposePartial( decomposed_partial_outer.First()));

      // no sub-label to extract
      if( decomposed_partial_inner.First().IsEmpty() || decomposed_partial_inner.Second().IsEmpty())
      {
        return LABEL;
      }

      // now perform the mapping
      storage::Vector< size_t> &partial_indices_outer( decomposed_partial_outer.Second());
      const storage::Vector< size_t> &partial_indices_inner( decomposed_partial_inner.Second());
      for
      (
        storage::Vector< size_t>::iterator
          itr_outer( partial_indices_outer.Begin()), itr_outer_end( partial_indices_outer.End());
        itr_outer != itr_outer_end;
        ++itr_outer
      )
      {
        BCL_Assert( *itr_outer < partial_indices_inner.GetSize(), "Invalid partial!");
        *itr_outer = partial_indices_inner( *itr_outer);
      }

      // create the serialization handler for the partials' indices
      util::OwnPtr< io::SerializationInterface> partial_indices_outer_handler
      (
        io::Serialization::GetAgent( &partial_indices_outer)
      );

      // create the new label with the remapped indices
      storage::Vector< util::ObjectDataLabel> new_arguments( 2);
      new_arguments( 0) = decomposed_partial_inner.First();
      new_arguments( 1) = partial_indices_outer_handler->GetLabel();
      new_arguments( 1).SetValue( s_indices);
      util::ObjectDataLabel new_label( s_partial, new_arguments);
      if( new_arguments( 0).GetValue() == s_partial)
      {
        return CollapsePartials( new_label);
      }
      return new_label;
    }

    //! @brief retrieve the indices of a sub-property
    //! @param LABEL property of the form Partial(x,indices(y1,y2,...,yn))
    //! @return y1,y2,...yn of the property x, if it could be located
    storage::Vector< size_t> FeatureLabelSet::GetSubPropertyIndices( const util::ObjectDataLabel &LABEL) const
    {
      // decompose the partial into partial and indices
      storage::Pair< util::ObjectDataLabel, storage::Vector< size_t> >
        decomposed_partial( DecomposePartial( LABEL));

      if( decomposed_partial.First().IsEmpty())
      {
        if( !m_ImplementationInterface.IsDefined())
        {
          BCL_Exit( LABEL.ToString() + " was not found in feature label set " + GetString(), -1);
        }
        util::Implementation< util::ImplementationInterface> impl_copy( m_ImplementationInterface.HardCopy());
        std::stringstream err_stream;
        BCL_Assert
        (
          impl_copy->TryRead( LABEL, err_stream),
          LABEL.ToString() + " was not found in feature label set " + GetString()
          + " and could not be standardized by " + impl_copy.GetString() + " because: "
          + err_stream.str()
        );
        util::ObjectDataLabel new_label( impl_copy->GetLabel());
        if( new_label != LABEL)
        {
          return GetPropertyIndices( new_label);
        }
      }

      // look for the label in the map
      storage::Map< util::ObjectDataLabel, math::Range< size_t> >::const_iterator itr_map
      (
        m_PropertiesToRanges.Find( decomposed_partial.First())
      );
      storage::Vector< size_t> indices_main_property;
      if( itr_map != m_PropertiesToRanges.End())
      {
        indices_main_property = GetPropertyIndices( decomposed_partial.First());

        // reorder to whatever order the indices were given in
        indices_main_property.Reorder( decomposed_partial.Second());
      }
      else
      {
        // this guarantees that any sub-features are located

        // create a feature label set for the requested feature, fully split
        FeatureLabelSet feature_for_label
        (
          "",
          storage::Vector< util::ObjectDataLabel>( size_t( 1), LABEL),
          storage::Vector< size_t>( size_t( 1), decomposed_partial.Second().GetSize()),
          m_ImplementationInterface
        );
        if( feature_for_label.GetSize() > size_t( 1))
        {
          feature_for_label = feature_for_label.SplitFeatureLabelSet();
        }

        if
        (
          math::Statistics::MaximumValue( m_PropertiesSizes.Begin(), m_PropertiesSizes.End()) > size_t( 1)
          || math::Statistics::MaximumValue( feature_for_label.GetPropertySizes().Begin(), feature_for_label.GetPropertySizes().End()) > size_t( 1)
        )
        {
          // create a split feature label set
          if( !m_SplitFeatures.IsDefined())
          {
            m_SplitFeatures = util::OwnPtr< FeatureLabelSet>( SplitFeatureLabelSet().Clone());
          }

          // get the common features
          for
          (
            storage::Vector< util::ObjectDataLabel>::const_iterator
              itr_features_label( feature_for_label.GetMemberLabels().Begin()),
              itr_features_label_end( feature_for_label.GetMemberLabels().End());
            itr_features_label != itr_features_label_end;
            ++itr_features_label
          )
          {
            indices_main_property.Append( m_SplitFeatures->GetPropertyIndices( *itr_features_label));
          }
        }
        else
        {
          BCL_Exit( "Could not locate property: " + LABEL.ToString() + " within " + GetString(), -1);
        }
      }

      return indices_main_property;
    }

  } // namespace model
} // namespace bcl
