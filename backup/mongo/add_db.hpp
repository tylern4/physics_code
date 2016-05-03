/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef WVSQ2_H_GUARD
#define WVSQ2_H_GUARD

#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>

void add_db(double W, double Q2) {
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{}};

    bsoncxx::builder::stream::document document{};

    auto collection = conn["physics"]["v2"];
    document << "W" << W;
    document << "Q2" << Q2;

    collection.insert_one(document.view());

}



#endif