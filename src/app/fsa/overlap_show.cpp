#include "overlap_show.hpp"

#include <cassert>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <sstream>
#include <iostream>

#include "utility.hpp"

OverlapShow::OverlapShow() {
   
}


bool OverlapShow::ParseArgument(int argc, const char *const argv[]) {
    ArgumentParser &&ap = GetArgumentParser();

    if (ap.ParseArgument(argc, argv)) {


       type_ = ol_store_.DetectFileType(ifname_);
        return true;
    }
    else {
        return false;
    }
}

ArgumentParser OverlapShow::GetArgumentParser() {
    ArgumentParser ap;

    ap.AddPositionOption(ifname_, "ifname", "overlap file path");
    ap.AddPositionOption(read_name_, "read_name", "read name");
    return ap;
}

void OverlapShow::Usage() {
    std::cout << GetArgumentParser().Usage();
}


void OverlapShow::Run() {
    LOG(INFO)("Load overlap file: %s", ifname_.c_str());
    LoadOverlaps(ifname_);

    //int id = type_ != "m4" ? ol_store_.GetReadStore().NameToId(read_name_) : std::atoi(read_name_.c_str());
    int id = 23;
    int max_overhang_ = 200;
    for (const auto& o : ol_store_.Get()) {
        if (o.Location(max_overhang_) != Overlap::Loc::Abnormal && (o.a_.id == id || o.b_.id == id)) {
            std::cout << ToString(o);
        }
    }
    LOG(INFO)("Load overlap file: %s", ifname_.c_str());
}

void OverlapShow::LoadOverlaps(const std::string &fname) {

    ol_store_.Load(ifname_);
    

}

std::string OverlapShow::ToString(const Overlap &o) const {
    auto t = ol_store_.DetectFileType(ifname_);
    if (t == "m4") return ol_store_.ToM4Line(o);
    else if (t == "m4a") return ol_store_.ToM4aLine(o);
    else if (t == "paf") return ol_store_.ToPafLine(o);
    else return ol_store_.ToM4Line(o);
}