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
#include "chemistry/bcl_chemistry_fragment_react.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // initialize static
    sched::Mutex &FragmentReact::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentReact::FragmentReact
    (
      const ReactionSearch &CATALOG
    ) :
      m_OutputMatchedReagents( false),
      m_ReactionWorker( ReactionWorker()),
      m_ReactionSearch( CATALOG)
    {
    }

    //! @brief reaction worker constructor
    FragmentReact::FragmentReact
    (
      const ReactionWorker &WORKER
    ) :
      m_OutputMatchedReagents( false),
      m_ReactionWorker( WORKER),
      m_ReactionSearch( ReactionSearch())
    {
    }

    //! @brief full constructor
    FragmentReact::FragmentReact
    (
      const ReactionWorker &WORKER,
      const ReactionSearch &CATALOG
    ) :
      m_OutputMatchedReagents( false),
      m_ReactionWorker( WORKER),
      m_ReactionSearch( CATALOG)
    {
    }

    //! @brief clone constructor
    FragmentReact *FragmentReact::Clone() const
    {
      return new FragmentReact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentReact::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentReact::GetAlias() const
    {
      static const std::string s_name( "React");
      return s_name;
    }

    //! @brief returns the reaction search object
    ReactionSearch &FragmentReact::GetReactionSearch() const
    {
      return m_ReactionSearch;
    }

    //! @brief returns the reaction worker object
    ReactionWorker &FragmentReact::GetReactionWorker() const
    {
      return m_ReactionWorker;
    }

    //! @brief returns the reactions
    const ReactionEnsemble &FragmentReact::GetReactions() const
    {
      return m_Reactions;
    }

    //! @brief returns the reagents
    const FragmentEnsemble &FragmentReact::GetReagents() const
    {
      return m_Reagents;
    }

    //! @brief returns the allowed reactant positions for the target molecule
    const storage::Vector< size_t> &FragmentReact::GetTargetReactantPositions() const
    {
      return m_TargetReactantPositions;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief set the reaction search object
    void FragmentReact::SetReactionSearch( const ReactionSearch &REACTION_SEARCH)
    {
      m_ReactionSearch = REACTION_SEARCH;
    }

    //! @brief set the reaction worker
    void FragmentReact::SetReactionWorker( const ReactionWorker &REACTION_WORKER)
    {
      m_ReactionWorker = REACTION_WORKER;
    }

    //! @brief set the reactions
    void FragmentReact::SetReactions( const ReactionEnsemble &REACTIONS)
    {
      m_Reactions = REACTIONS;
    }

    //! @brief set the reagents
    void FragmentReact::SetReagents( const FragmentEnsemble &REAGENTS)
    {
      m_Reagents = REAGENTS;
    }

    //! @brief set the allowed reactant positions for the target molecule
    void FragmentReact::SetTargetReactantPositions( const storage::Vector< size_t> &TARGET_REACTANT_POS)
    {
      m_TargetReactantPositions = TARGET_REACTANT_POS;
    }

    //! @brief react a molecule at a random reactive center
    storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> FragmentReact::ReactRandom
    ( 
      const FragmentComplete &MOL
    ) const
    {
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res;

      // This SiPtrVector< const FragmentComplete> corresponds to the different reactant positions for a single reaction
      // (i.e. as opposed to all the reactants available for one position of a reaction)
      // so looping over it is equivalent to looping over the positions of reactants in the reaction
      storage::Pair< util::SiPtr< const ReactionComplete>, util::SiPtrVector< const FragmentComplete> >
        random_rxn_rxts( m_ReactionSearch.ChooseRandomRxnAndReactants( MOL));
      storage::Vector< util::SiPtrVector< const FragmentComplete> > rxn_search( m_ReactionSearch.GetAvailableReactants( random_rxn_rxts.First()));

      if( !random_rxn_rxts.First().IsDefined())
      {
        BCL_MessageVrb( "ReactRandom: Could not find a suitable reaction for molecule " + RemoveWhitespace( MOL.GetName()));
        return res;
      }

      storage::Vector< size_t> mol_rxt_pos( m_ReactionWorker.MatchesReactants( MOL, *random_rxn_rxts.First()));

      if( mol_rxt_pos.IsEmpty())
      {
        BCL_MessageVrb( "ReactRandom: Mol did not actually match any reactants");
        return res;
      }
      mol_rxt_pos.Shuffle();
      size_t &picked( mol_rxt_pos( 0));

      FragmentEnsemble reactants;
      for( size_t i( 0), l( random_rxn_rxts.Second().GetSize()); i < l; ++i)
      {
        if( i == picked)
        {
          reactants.PushBack( MOL);
        }
        else
        {
          if( !random_rxn_rxts.Second()( i).IsDefined())
          {
            return res;
          }
          reactants.PushBack( *random_rxn_rxts.Second()( i));
        }
      }

      // get position of primary reagent
      size_t n_pos( rxn_search.GetSize());
      size_t pos( 0);
      for( ; pos < n_pos; ++pos)
      {
        // if the main reagent matches multiple positions, give preference to earlier position
        if( m_ReactionWorker.MatchesReactantIndex( MOL, *random_rxn_rxts.First(), pos, m_TargetReactantPositions))
        {
          break;
        }
      }

      res.First() = random_rxn_rxts.First();
      res.Second() = m_ReactionWorker.ExecuteReaction( *random_rxn_rxts.First(), reactants, MOL, pos);
      return res;
    }

    //! @brief react a molecule at all reactive centers with all reagents for a single specified reaction
    storage::Pair
    <
      storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> >,
      storage::Vector< storage::Vector< std::string> >
    >
    FragmentReact::ReactExhaustiveOneReaction
    (
      const FragmentComplete &MOL,
      const ReactionComplete &RXN
    ) const
    {
      storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> > res_v;
      storage::Vector< storage::Vector< std::string> > rxn_info_v;

      // for each available reaction, find the available reactants
      //      m_ReactionSearch.Initialize();
      util::SiPtrVector< const ReactionComplete> avail_reactions( m_ReactionSearch.FindReactions( MOL));
      for
      (
          auto itr( avail_reactions.Begin()), itr_end( avail_reactions.End());
          itr != itr_end;
          ++itr
      )
      {
        if( ( *itr)->GetDescription() == RXN.GetDescription())
        {
          BCL_MessageStd( "CURRENT RXN: " + util::Format()( RXN.GetDescription()));
          if
          (
              !m_ReactionSearch.GetAvailableReactants( *itr).IsEmpty() &&
              ( *itr)->GetNumberReactants() == m_ReactionSearch.GetAvailableReactants( *itr).GetSize()
          )
          {
            size_t n_rxt_pos( m_ReactionSearch.GetAvailableReactants( *itr).GetSize());
            BCL_MessageStd( "Number of reactant positions in reaction: " + util::Format()( n_rxt_pos));

            // TODO: refactor to call recursively and support arbitrary number of reagents
            // call helper for appropriate number of reactant positions
            // one reactant position
            if( n_rxt_pos == size_t( 1))
            {
              auto prod( ReactExhaustiveOneReactantHelper( *itr, MOL));
              res_v.PushBack( prod.First());
              rxn_info_v.PushBack( prod.Second());
            }
            // two reactant positions
            else if( n_rxt_pos == size_t( 2))
            {
              auto prod( ReactExhaustiveTwoReactantHelper( *itr, MOL));
              res_v.PushBack( prod.First());
              rxn_info_v.PushBack( prod.Second());
            }
            // three reactant positions
            else if( n_rxt_pos == size_t( 3))
            {
              auto prod( ReactExhaustiveThreeReactantHelper( *itr, MOL));
              res_v.PushBack( prod.First());
              rxn_info_v.PushBack( prod.Second());
            }
            // four reactant positions
            else if( n_rxt_pos == size_t( 4))
            {
              auto prod( ReactExhaustiveFourReactantHelper( *itr, MOL));
              res_v.PushBack( prod.First());
              rxn_info_v.PushBack( prod.Second());
            }
          }
        }
      }
      return std::make_pair( res_v, rxn_info_v);
    }

    //! @brief react a molecule at all reactive centers with all reagents
    storage::Pair
    <
      storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> >,
      storage::Vector< storage::Vector< std::string> >
    >
    FragmentReact::ReactExhaustive
    (
      const FragmentComplete &MOL
    ) const
    {
      storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> > res_v;
      storage::Vector< storage::Vector< std::string> > rxn_info_v;

      // find available reactions
      util::SiPtrVector< const ReactionComplete> avail_reactions( m_ReactionSearch.FindReactions( MOL));

      // if there are no available reactions, then bail
      if( avail_reactions.IsEmpty())
      {
        BCL_MessageStd( "ReactExhaustive: Could not find any suitable reactions for molecule " + RemoveWhitespace( MOL.GetName()));
        return std::make_pair( res_v, rxn_info_v);
      }

      // go over each available reaction
      for
      (
          auto rxns_itr( avail_reactions.Begin()), rxns_itr_end( avail_reactions.End());
          rxns_itr != rxns_itr_end;
          ++rxns_itr
      )
      {
        // for each available reaction, find the available reactants
        BCL_MessageStd( "CURRENT RXN: " + util::Format()( ( *rxns_itr)->GetDescription()));
        if
        (
            !m_ReactionSearch.GetAvailableReactants( *rxns_itr).IsEmpty() &&
            ( *rxns_itr)->GetNumberReactants() == m_ReactionSearch.GetAvailableReactants( *rxns_itr).GetSize()
        )
        {
          size_t n_rxt_pos( m_ReactionSearch.GetAvailableReactants( *rxns_itr).GetSize());
          BCL_MessageStd( "Number of reactant positions in reaction: " + util::Format()( n_rxt_pos));

          // TODO: refactor to make recursive and generalized to N reactants
          // call helper for apprioriate number of reactant positions
          // one reactant position
          if( n_rxt_pos == size_t( 1))
          {
            auto prod( ReactExhaustiveOneReactantHelper( *rxns_itr, MOL));
            res_v.PushBack( prod.First());
            rxn_info_v.PushBack( prod.Second());
          }
          // two reactant positions
          else if( n_rxt_pos == size_t( 2))
          {
            auto prod( ReactExhaustiveTwoReactantHelper( *rxns_itr, MOL));
            res_v.PushBack( prod.First());
            rxn_info_v.PushBack( prod.Second());
          }
          // three or more reactant positions
          else if( n_rxt_pos == size_t( 3))
          {
            auto prod( ReactExhaustiveThreeReactantHelper( *rxns_itr, MOL));
            res_v.PushBack( prod.First());
            rxn_info_v.PushBack( prod.Second());
          }
          else if( n_rxt_pos == size_t( 4))
          {
            auto prod( ReactExhaustiveFourReactantHelper( *rxns_itr, MOL));
            res_v.PushBack( prod.First());
            rxn_info_v.PushBack( prod.Second());
          }
          // we currently do not support more than 4 reagents per reaction
          else
          {
            continue;
          }
        }
      }
      return std::make_pair( res_v, rxn_info_v);
    }

    //! @brief worker function for three reactant reactions
    //! @param RXN the reaction to be performed
    //! @param MOL the molecule on which to perform the reaction
    //! @return returns a pair where the first element is a pair
    //! containing the full ReactionComplete and the ensemble of products,
    //! and the second element is a vector of strings corresponding to
    //! each product in the output ensemble indicating the reaction and reactants
    //! used to produce said product
    storage::Pair
    <
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
      storage::Vector< std::string>
    >
    FragmentReact::ReactExhaustiveOneReactantHelper
    (
      const util::SiPtr< const ReactionComplete> &RXN,
      const FragmentComplete &MOL
    ) const
    {
      // initialize output
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res;
      storage::Vector< std::string> rxn_info;

      // initialize reaction/reactants pair
      storage::Pair< util::SiPtr< const ReactionComplete>, util::SiPtrVector< const FragmentComplete> > rxn_rxt_pair;

      // search reaction to see if MOL could be a reactant
      storage::Vector< size_t> mol_rxt_pos( m_ReactionWorker.MatchesReactants( MOL, *RXN, m_TargetReactantPositions));

      // make sure we have some matches, else continue to next reaction/reactant pair
      if( mol_rxt_pos.IsEmpty())
      {
        BCL_MessageStd
        (
          "ReactExhaustive: MOL is not a reactant for reaction " + util::Format()( *RXN)
        );
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }

      // this should be a vector of size 1 in this function corresponding to a single reactant position
      storage::Vector< util::SiPtrVector< const FragmentComplete>> rxn_search( m_ReactionSearch.GetAvailableReactants( RXN));
      BCL_MessageStd( "Expected 1 reactant position(s) and found " + util::Format()( rxn_search.GetSize()) + " reactant position(s)");
      for
      (
          auto rxt_pos_itr( rxn_search.Begin()),
          rxt_itr_end( rxn_search.End());
          rxt_pos_itr != rxt_itr_end;
          ++rxt_pos_itr
      )
      {
        // reaction and vector of reactants for this reactant position
        rxn_rxt_pair = std::make_pair( RXN, util::SiPtrVector< const FragmentComplete>( *rxt_pos_itr));

        // make sure there is a reactant for the reaction
        if( rxn_rxt_pair.Second().GetSize())
        {
          FragmentEnsemble reactants;
          reactants.PushBack( MOL);

          // store rxn and reaction products
          res.First() = RXN;
          FragmentEnsemble products( m_ReactionWorker.ExecuteReaction( *RXN, reactants, MOL, size_t( 0)));

          // save products
          for
          (
              auto prod_itr( products.Begin()), prod_itr_end( products.End());
              prod_itr != prod_itr_end;
              ++prod_itr
          )
          {
            // make reactant string label
            std::string rxt_label;
            size_t rxt_index( 0);
            for
            (
                auto rxt_itr( reactants.Begin()), rxt_itr_end( reactants.End());
                rxt_itr != rxt_itr_end;
                ++rxt_itr, ++rxt_index
            )
            {
              std::istringstream rxt_ss( rxt_itr->GetName());
              std::string rxt_line;
              std::getline( rxt_ss, rxt_line);
              rxt_label.append( RemoveWhitespace( rxt_line));
              if( rxt_index < reactants.GetSize() - 1)
              {
                rxt_label.append( ",");
              }
            }

            // store reaction as property
            prod_itr->StoreProperty
            (
              "Reaction",
              util::Format()( RemoveWhitespace( RXN->GetDescription())) +
              "," + util::Format()( rxt_label)
            );
            rxn_info.PushBack( prod_itr->GetMDLProperty( "Reaction"));
            res.Second().PushBack( *prod_itr);
          }
        }
        // no reactants for this position
        else
        {
          return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
        }
      }
      BCL_MessageStd( "Generated " + util::Format()( res.Second().GetSize()) + " product(s) with this reaction!");
      if( res.Second().GetSize())
      {
        return std::make_pair
            (
              storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>( std::make_pair( res.First(), res.Second())),
              rxn_info
            );
      }
      else
      {
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }
    }

    //! @brief worker function for three reactant reactions
    //! @param RXN the reaction to be performed
    //! @param MOL the molecule on which to perform the reaction
    //! @return returns a pair where the first element is a pair
    //! containing the full ReactionComplete and the ensemble of products,
    //! and the second element is a vector of strings corresponding to
    //! each product in the output ensemble indicating the reaction and reactants
    //! used to produce said product
    storage::Pair
    <
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
      storage::Vector< std::string>
    >
    FragmentReact::ReactExhaustiveTwoReactantHelper
    (
      const util::SiPtr< const ReactionComplete> &RXN,
      const FragmentComplete &MOL
    ) const
    {
      // debug output
      io::OFStream debug_mdl_out, rxn_pos_0, rxn_pos_1;
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::MustOpenOFStream( debug_mdl_out, "ReactExhaustive2.sdf");
        io::File::MustOpenOFStream( rxn_pos_0, "rxn_pos_0.sdf");
        io::File::MustOpenOFStream( rxn_pos_1, "rxn_pos_1.sdf");
        GetMutex().Unlock();
      }

      // initialize output
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res;
      storage::Vector< std::string> rxn_info;

      // initialize reaction/reactants pair
      storage::Pair< util::SiPtr< const ReactionComplete>, util::SiPtrVector< const FragmentComplete> > rxn_rxt_pair;

      // search reaction to see if MOL could be a reactant
      // this tells us which positions in the reaction the molecule could act as a reactant
      storage::Vector< size_t> mol_rxt_pos( m_ReactionWorker.MatchesReactants( MOL, *RXN, m_TargetReactantPositions));

      // make sure we have some matches, else continue to next reaction/reactant pair
      if( mol_rxt_pos.IsEmpty())
      {
        BCL_MessageStd
        (
          "ReactExhaustive: MOL is not a reactant for reaction " + util::Format()( *RXN)
        );
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }

      // this should be a vector of size 2 in this function corresponding to 2 reactant positions
      storage::Vector< util::SiPtrVector< const FragmentComplete>> rxn_search( m_ReactionSearch.GetAvailableReactants( RXN));

      // debug
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        for( auto blah_0( rxn_search( 0).Begin()), blah_0_end( rxn_search( 0).End()); blah_0 != blah_0_end; ++blah_0)
        {
          ( *blah_0)->WriteMDL( rxn_pos_0);
        }
        for( auto blah_1( rxn_search( 1).Begin()), blah_1_end( rxn_search( 1).End()); blah_1 != blah_1_end; ++blah_1)
        {
          ( *blah_1)->WriteMDL( rxn_pos_1);
        }
        GetMutex().Unlock();
      }

      size_t n_pos( rxn_search.GetSize());
      BCL_MessageStd( "Expected 2 reactant position(s) and found " + util::Format()( n_pos) + " reactant position(s)");

      // fill out our reactant pairs for this reaction
      storage::Vector< FragmentEnsemble> all_products;

      // reactant positions
      // fill each reactant position we can with MOL and enumerate the paired reactants
      for( size_t pos( 0); pos < n_pos; ++pos)
      {
        // check to see if MOL can be a reactant in this position
        if( mol_rxt_pos.Find( pos) >= mol_rxt_pos.GetSize())
        {
          continue;
        }
        // enumerate the reactants that can go in the paired position
        for( size_t other_pos( 0); other_pos < n_pos; ++other_pos)
        {
          // cannot simultaneously occupy same reactant position
          if( other_pos == pos)
          {
            continue;
          }
          // set up our vector indexing reactants by their position in the reaction
          storage::Vector< FragmentComplete> partial_rxts( n_pos);

          // assign the position of our molecule of interest
          partial_rxts( pos) = MOL;

          // now enumerate the reactants in the original SiPtrVector< FragmentComplete> at the correct position
          if( rxn_search( other_pos).GetSize())
          {
            for
            (
                auto rxn_search_itr( rxn_search( other_pos).Begin()), rxn_search_itr_end( rxn_search( other_pos).End());
                rxn_search_itr != rxn_search_itr_end;
                ++rxn_search_itr
            )
            {
              // copy our partial reactants vector
              storage::Vector< FragmentComplete> rxts( partial_rxts);
              // add the thing
              rxts( other_pos) = **rxn_search_itr;

              // make our ensemble for reactants
              FragmentEnsemble reactants( storage::List< FragmentComplete>( rxts.Begin(), rxts.End()));

              // debug write out reactants
              if( m_OutputMatchedReagents)
              {
                GetMutex().Lock();
                for( auto itr( reactants.Begin()); itr != reactants.End(); ++itr)
                {
                  itr->WriteMDL( debug_mdl_out);
                }
                GetMutex().Unlock();
              }

              // perform reaction with this set of reactants
              FragmentEnsemble products( m_ReactionWorker.ExecuteReaction( *RXN, reactants, MOL, pos));

              // save reaction information
              for
              (
                  auto prod_itr( products.Begin()), prod_itr_end( products.End());
                  prod_itr != prod_itr_end;
                  ++prod_itr
              )
              {
                // make reactant string label
                std::string rxt_label;
                size_t rxt_index( 0);
                for
                (
                    auto rxt_itr( reactants.Begin()), rxt_itr_end( reactants.End());
                    rxt_itr != rxt_itr_end;
                    ++rxt_itr, ++rxt_index
                )
                {
                  std::istringstream rxt_ss( rxt_itr->GetName());
                  std::string rxt_line;
                  std::getline( rxt_ss, rxt_line);
                  rxt_label.append( RemoveWhitespace( rxt_line));
                  if( rxt_index < reactants.GetSize() - 1)
                  {
                    rxt_label.append( ",");
                  }
                }

                // store reaction as property
                prod_itr->StoreProperty
                (
                  "Reaction",
                  util::Format()( RemoveWhitespace( RXN->GetDescription())) +
                  "," + util::Format()( rxt_label)
                );
                rxn_info.PushBack( prod_itr->GetMDLProperty( "Reaction"));
              }
              all_products.PushBack( products);
            }
          }
        }
      }

      // close debug write
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::CloseClearFStream( debug_mdl_out);
        io::File::CloseClearFStream( rxn_pos_0);
        io::File::CloseClearFStream( rxn_pos_1);
        GetMutex().Unlock();
      }

      // finalize our output
      res.First() = RXN;
      for( size_t p( 0); p < all_products.GetSize(); ++p)
      {
        for
        (
            auto prod_itr( all_products( p).Begin()), prod_itr_end( all_products( p).End());
            prod_itr != prod_itr_end;
            ++prod_itr
        )
        {
          res.Second().PushBack( *prod_itr);
        }
      }
      BCL_MessageStd( "Generated " + util::Format()( res.Second().GetSize()) + " product(s) with this reaction!");
      if( res.Second().GetSize())
      {
        return std::make_pair
            (
              storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>( std::make_pair( res.First(), res.Second())),
              rxn_info
            );
      }
      else
      {
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }
    }

    //! @brief worker function for three reactant reactions
    //! @param RXN the reaction to be performed
    //! @param MOL the molecule on which to perform the reaction
    //! @return returns a pair where the first element is a pair
    //! containing the full ReactionComplete and the ensemble of products,
    //! and the second element is a vector of strings corresponding to
    //! each product in the output ensemble indicating the reaction and reactants
    //! used to produce said product
    storage::Pair
    <
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
      storage::Vector< std::string>
    >
    FragmentReact::ReactExhaustiveThreeReactantHelper
    (
      const util::SiPtr< const ReactionComplete> &RXN,
      const FragmentComplete &MOL
    ) const
    {
      // debug output
      io::OFStream debug_mdl_out, rxn_pos_0, rxn_pos_1, rxn_pos_2;
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::MustOpenOFStream( debug_mdl_out, "ReactExhaustive3.sdf");
        io::File::MustOpenOFStream( rxn_pos_0, "rxn_pos_0.sdf");
        io::File::MustOpenOFStream( rxn_pos_1, "rxn_pos_1.sdf");
        io::File::MustOpenOFStream( rxn_pos_2, "rxn_pos_2.sdf");
        GetMutex().Unlock();
      }

      // initialize output
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res;
      storage::Vector< std::string> rxn_info;

      // initialize reaction/reactants pair
      storage::Pair< util::SiPtr< const ReactionComplete>, util::SiPtrVector< const FragmentComplete> > rxn_rxt_pair;

      // search reaction to see if MOL could be a reactant
      // this tells us which positions in the reaction the molecule could act as a reactant
      storage::Vector< size_t> mol_rxt_pos( m_ReactionWorker.MatchesReactants( MOL, *RXN, m_TargetReactantPositions));

      // make sure we have some matches, else continue to next reaction/reactant pair
      if( mol_rxt_pos.IsEmpty())
      {
        BCL_MessageStd
        (
          "ReactExhaustive: MOL is not a reactant for reaction " + util::Format()( *RXN)
        );
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }

      // this should be a vector of size 3 in this function corresponding to 3 reactant positions
      storage::Vector< util::SiPtrVector< const FragmentComplete>> rxn_search( m_ReactionSearch.GetAvailableReactants( RXN));

      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        for( auto blah_0( rxn_search( 0).Begin()), blah_0_end( rxn_search( 0).End()); blah_0 != blah_0_end; ++blah_0)
        {
          ( *blah_0)->WriteMDL( rxn_pos_0);
        }
        for( auto blah_1( rxn_search( 1).Begin()), blah_1_end( rxn_search( 1).End()); blah_1 != blah_1_end; ++blah_1)
        {
          ( *blah_1)->WriteMDL( rxn_pos_1);
        }
        for( auto blah_2( rxn_search( 2).Begin()), blah_2_end( rxn_search( 2).End()); blah_2 != blah_2_end; ++blah_2)
        {
          ( *blah_2)->WriteMDL( rxn_pos_2);
        }
        GetMutex().Unlock();
      }

      size_t n_pos( rxn_search.GetSize());
      BCL_MessageStd( "Expected 3 reactant position(s) and found " + util::Format()( n_pos) + " reactant position(s)");

      // fill out our reactant pairs for this reaction
      storage::Vector< FragmentEnsemble> all_products;

      // iterate over molecules in first position
      for
      (
          auto rxn_search_itr( rxn_search( 0).Begin()), rxn_search_itr_end( rxn_search( 0).End());
          rxn_search_itr != rxn_search_itr_end;
          ++rxn_search_itr
      )
      {
        // iterate over molecules in second position
        for
        (
            auto rxn_search_itr_b( rxn_search( 1).Begin()), rxn_search_itr_b_end( rxn_search( 1).End());
            rxn_search_itr_b != rxn_search_itr_b_end;
            ++rxn_search_itr_b
        )
        {
          // iterate over molecules in third position
          for
          (
              auto rxn_search_itr_c( rxn_search( 2).Begin()), rxn_search_itr_c_end( rxn_search( 2).End());
              rxn_search_itr_c != rxn_search_itr_c_end;
              ++rxn_search_itr_c
          )
          {
            // make our ensemble for reactants
            storage::List< FragmentComplete> reactant_list;
            reactant_list.InsertElement( **rxn_search_itr);
            reactant_list.InsertElement( **rxn_search_itr_b);
            reactant_list.InsertElement( **rxn_search_itr_c);
            FragmentEnsemble reactants( reactant_list);

            // debug write out reactants
            if( m_OutputMatchedReagents)
            {
              GetMutex().Lock();
              for( auto itr( reactants.Begin()); itr != reactants.End(); ++itr)
              {
                itr->WriteMDL( debug_mdl_out);
              }
              GetMutex().Unlock();
            }

            // get position of primary reagent
            size_t pos( 0);
            for( ; pos < n_pos; ++pos)
            {
              // if the main reagent matches multiple positions, give preference to earlier position
              if( m_ReactionWorker.MatchesReactantIndex( MOL, *RXN, pos, m_TargetReactantPositions))
              {
                break;
              }
            }
            FragmentEnsemble products( m_ReactionWorker.ExecuteReaction( *RXN, reactants, MOL, pos));

            // save reaction information
            for
            (
                auto prod_itr( products.Begin()), prod_itr_end( products.End());
                prod_itr != prod_itr_end;
                ++prod_itr
            )
            {
              // make reactant string label
              std::string rxt_label;
              size_t rxt_index( 0);
              for
              (
                  auto rxt_itr( reactants.Begin()), rxt_itr_end( reactants.End());
                  rxt_itr != rxt_itr_end;
                  ++rxt_itr, ++rxt_index
              )
              {
                std::istringstream rxt_ss( rxt_itr->GetName());
                std::string rxt_line;
                std::getline( rxt_ss, rxt_line);
                rxt_label.append( RemoveWhitespace( rxt_line));
                if( rxt_index < reactants.GetSize() - 1)
                {
                  rxt_label.append( ",");
                }
              }

              // store reaction as property
              prod_itr->StoreProperty
              (
                "Reaction",
                util::Format()( RemoveWhitespace( RXN->GetDescription())) +
                "," + util::Format()( rxt_label)
              );
              rxn_info.PushBack( prod_itr->GetMDLProperty( "Reaction"));
            }
            all_products.PushBack( products);
          }
        }
      }

      // close debug write
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::CloseClearFStream( debug_mdl_out);
        io::File::CloseClearFStream( rxn_pos_0);
        io::File::CloseClearFStream( rxn_pos_1);
        io::File::CloseClearFStream( rxn_pos_2);
        GetMutex().Unlock();
      }

      // finalize our output
      res.First() = RXN;
      for( size_t p( 0); p < all_products.GetSize(); ++p)
      {
        for
        (
            auto prod_itr( all_products( p).Begin()), prod_itr_end( all_products( p).End());
            prod_itr != prod_itr_end;
            ++prod_itr
        )
        {
          res.Second().PushBack( *prod_itr);
        }
      }
      BCL_MessageStd( "Generated " + util::Format()( res.Second().GetSize()) + " product(s) with this reaction!");
      if( res.Second().GetSize())
      {
        return std::make_pair
            (
              storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>( std::make_pair( res.First(), res.Second())),
              rxn_info
            );
      }
      else
      {
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }
    }

    //! @brief worker function for three reactant reactions
    //! @param RXN the reaction to be performed
    //! @param MOL the molecule on which to perform the reaction
    //! @return returns a pair where the first element is a pair
    //! containing the full ReactionComplete and the ensemble of products,
    //! and the second element is a vector of strings corresponding to
    //! each product in the output ensemble indicating the reaction and reactants
    //! used to produce said product
    storage::Pair
    <
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
      storage::Vector< std::string>
    >
    FragmentReact::ReactExhaustiveFourReactantHelper
    (
      const util::SiPtr< const ReactionComplete> &RXN,
      const FragmentComplete &MOL
    ) const
    {
      // debug output
      io::OFStream debug_mdl_out, rxn_pos_0, rxn_pos_1, rxn_pos_2, rxn_pos_3;
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::MustOpenOFStream( debug_mdl_out, "ReactExhaustive3.sdf");
        io::File::MustOpenOFStream( rxn_pos_0, "rxn_pos_0.sdf");
        io::File::MustOpenOFStream( rxn_pos_1, "rxn_pos_1.sdf");
        io::File::MustOpenOFStream( rxn_pos_2, "rxn_pos_2.sdf");
        io::File::MustOpenOFStream( rxn_pos_3, "rxn_pos_3.sdf");
        GetMutex().Unlock();
      }

      // initialize output
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> res;
      storage::Vector< std::string> rxn_info;

      // initialize reaction/reactants pair
      storage::Pair< util::SiPtr< const ReactionComplete>, util::SiPtrVector< const FragmentComplete> > rxn_rxt_pair;

      // search reaction to see if MOL could be a reactant
      // this tells us which positions in the reaction the molecule could act as a reactant
      storage::Vector< size_t> mol_rxt_pos( m_ReactionWorker.MatchesReactants( MOL, *RXN, m_TargetReactantPositions));

      // make sure we have some matches, else continue to next reaction/reactant pair
      if( mol_rxt_pos.IsEmpty())
      {
        BCL_MessageStd
        (
          "ReactExhaustive: MOL is not a reactant for reaction " + util::Format()( *RXN)
        );
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }

      // this should be a vector of size 4 in this function corresponding to 4 reactant positions
      storage::Vector< util::SiPtrVector< const FragmentComplete>> rxn_search( m_ReactionSearch.GetAvailableReactants( RXN));

      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        for( auto blah_0( rxn_search( 0).Begin()), blah_0_end( rxn_search( 0).End()); blah_0 != blah_0_end; ++blah_0)
        {
          ( *blah_0)->WriteMDL( rxn_pos_0);
        }
        for( auto blah_1( rxn_search( 1).Begin()), blah_1_end( rxn_search( 1).End()); blah_1 != blah_1_end; ++blah_1)
        {
          ( *blah_1)->WriteMDL( rxn_pos_1);
        }
        for( auto blah_2( rxn_search( 2).Begin()), blah_2_end( rxn_search( 2).End()); blah_2 != blah_2_end; ++blah_2)
        {
          ( *blah_2)->WriteMDL( rxn_pos_2);
        }
        for( auto blah_3( rxn_search( 3).Begin()), blah_3_end( rxn_search( 3).End()); blah_3 != blah_3_end; ++blah_3)
        {
          ( *blah_3)->WriteMDL( rxn_pos_3);
        }
        GetMutex().Unlock();
      }

      size_t n_pos( rxn_search.GetSize());
      BCL_MessageStd( "Expected 4 reactant position(s) and found " + util::Format()( n_pos) + " reactant position(s)");

      // fill out our reactant pairs for this reaction
      storage::Vector< FragmentEnsemble> all_products;

      // iterate over molecules in first position
      for
      (
          auto rxn_search_itr( rxn_search( 0).Begin()), rxn_search_itr_end( rxn_search( 0).End());
          rxn_search_itr != rxn_search_itr_end;
          ++rxn_search_itr
      )
      {
        // iterate over molecules in second position
        for
        (
            auto rxn_search_itr_b( rxn_search( 1).Begin()), rxn_search_itr_b_end( rxn_search( 1).End());
            rxn_search_itr_b != rxn_search_itr_b_end;
            ++rxn_search_itr_b
        )
        {
          // iterate over molecules in third position
          for
          (
              auto rxn_search_itr_c( rxn_search( 2).Begin()), rxn_search_itr_c_end( rxn_search( 2).End());
              rxn_search_itr_c != rxn_search_itr_c_end;
              ++rxn_search_itr_c
          )
          {
            // iterate over molecules in fourth position
            for
            (
                auto rxn_search_itr_d( rxn_search( 3).Begin()), rxn_search_itr_d_end( rxn_search( 3).End());
                rxn_search_itr_d != rxn_search_itr_d_end;
                ++rxn_search_itr_d
            )
            {
              // make our ensemble for reactants
              storage::List< FragmentComplete> reactant_list;
              reactant_list.InsertElement( **rxn_search_itr);
              reactant_list.InsertElement( **rxn_search_itr_b);
              reactant_list.InsertElement( **rxn_search_itr_c);
              reactant_list.InsertElement( **rxn_search_itr_d);
              FragmentEnsemble reactants( reactant_list);

              // debug write out reactants
              if( m_OutputMatchedReagents)
              {
                GetMutex().Lock();
                for( auto itr( reactants.Begin()); itr != reactants.End(); ++itr)
                {
                  itr->WriteMDL( debug_mdl_out);
                }
                GetMutex().Unlock();
              }

              // get position of primary reagent
              size_t pos( 0);
              for( ; pos < n_pos; ++pos)
              {
                // if the main reagent matches multiple positions, give preference to earlier position
                if( m_ReactionWorker.MatchesReactantIndex( MOL, *RXN, pos, m_TargetReactantPositions))
                {
                  break;
                }
              }
              FragmentEnsemble products( m_ReactionWorker.ExecuteReaction( *RXN, reactants, MOL, pos));

              // save reaction information
              for
              (
                  auto prod_itr( products.Begin()), prod_itr_end( products.End());
                  prod_itr != prod_itr_end;
                  ++prod_itr
              )
              {
                // make reactant string label
                std::string rxt_label;
                size_t rxt_index( 0);
                for
                (
                    auto rxt_itr( reactants.Begin()), rxt_itr_end( reactants.End());
                    rxt_itr != rxt_itr_end;
                    ++rxt_itr, ++rxt_index
                )
                {
                  std::istringstream rxt_ss( rxt_itr->GetName());
                  std::string rxt_line;
                  std::getline( rxt_ss, rxt_line);
                  rxt_label.append( RemoveWhitespace( rxt_line));
                  if( rxt_index < reactants.GetSize() - 1)
                  {
                    rxt_label.append( ",");
                  }
                }

                // store reaction as property
                prod_itr->StoreProperty
                (
                  "Reaction",
                  util::Format()( RemoveWhitespace( RXN->GetDescription())) +
                  "," + util::Format()( rxt_label)
                );
                rxn_info.PushBack( prod_itr->GetMDLProperty( "Reaction"));
              }
              all_products.PushBack( products);
            }
          }
        }
      }

      // close debug write
      if( m_OutputMatchedReagents)
      {
        GetMutex().Lock();
        io::File::CloseClearFStream( debug_mdl_out);
        io::File::CloseClearFStream( rxn_pos_0);
        io::File::CloseClearFStream( rxn_pos_1);
        io::File::CloseClearFStream( rxn_pos_2);
        io::File::CloseClearFStream( rxn_pos_3);
        GetMutex().Unlock();
      }

      // finalize our output
      res.First() = RXN;
      for( size_t p( 0); p < all_products.GetSize(); ++p)
      {
        for
        (
            auto prod_itr( all_products( p).Begin()), prod_itr_end( all_products( p).End());
            prod_itr != prod_itr_end;
            ++prod_itr
        )
        {
          res.Second().PushBack( *prod_itr);
        }
      }
      BCL_MessageStd( "Generated " + util::Format()( res.Second().GetSize()) + " product(s) with this reaction!");
      if( res.Second().GetSize())
      {
        return std::make_pair
            (
              storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>( std::make_pair( res.First(), res.Second())),
              rxn_info
            );
      }
      else
      {
        return storage::Pair< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>, storage::Vector< std::string>>();
      }
    }

    //! @brief react a molecule with a center containing a specified set of atoms
    //! @param REQUIRE_ALL whether all atom indices must be in a reactive center (true), or
    //!        if only a single atom must be in a reactive center (false)
    FragmentEnsemble FragmentReact::ReactCenterContaining
    ( 
      const FragmentComplete &MOL, 
      storage::Vector< size_t> &ATOM_INDICES,
      const bool &REQUIRE_ALL
    ) const
    {
      return FragmentEnsemble();
    }

    //! @brief remove whitespace (via isspace) from a string
     //! @param STR the string to remove whitespace from
     //! @return STR without any whitespace
     std::string FragmentReact::RemoveWhitespace( const std::string &STR) const
     {
       std::string str;
       str.reserve( STR.length());
       for
       (
         std::string::const_iterator itr( STR.begin()), itr_end( STR.end());
         itr != itr_end;
         ++itr
       )
       {
         if( !isspace( *itr))
         {
           str.push_back( *itr);
         }
       }
       return str;
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentReact::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentReact::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentReact::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "React fragments to generate products according to a RXN file and input reagents"
      );

      parameters.AddInitializer
      (
        "output_matched_reagents",
        "in addition to standard output, write SDFs corresponding to each reactant position that "
        "contain the input reagents matched to that position; this is useful for debugging RXN files; "
        "no effect on single component reactions or random reactions",
        io::Serialization::GetAgent( &m_OutputMatchedReagents),
        "false"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentReact::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
