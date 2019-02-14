#ifndef FSA_READ_STORE_HPP
#define FSA_READ_STORE_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include "sequence.hpp"

class ReadStore {
public:
    struct Item {
        Item(){}
        Item(const std::string &s, SeqReader::ItemId i, SeqReader* r) : seq(s), id(i), reader(r) {}
        std::string seq;
        SeqReader::ItemId id {-1};
        SeqReader* reader {nullptr};
    };

    void SetNameToId(const std::string &name, Seq::Id id);
    int NameToId(const std::string &name);
    std::string IdToName(Seq::Id id);
    const std::string& GetSeq(Seq::Id id) const { LoadItem(items_[id]); return items_[id].seq;  }
    const std::string& GetSeq(const std::string &name) { return GetSeq(names_to_ids_[name]); }

    std::string GetSeq(const Seq::Area& sa);
    size_t GetSeqLength(Seq::Id id) const { return GetSeq(id).size(); }

    void SaveIdToName(const std::string& fname) const;
    std::array<Seq::Id, 2> GetIdRange() const { return std::array<Seq::Id, 2>{0, (int)names_.size()}; }

    void Load(const std::string &fname, const std::string &type="", int mode=0);
    void LoadFasta(const std::string &fname, int mode=0);
    void LoadFastq(const std::string &fname, int mode=0);
    void LoadFofn(const std::string &fname, int mode=0);
    void LoadTxt(const std::string &fname, int mode=0) { LoadFofn(fname, mode); }
    void LoadItem(Item &item) const {
        if (item.seq.empty() && item.reader != nullptr) {
            SeqReader::Item i;
            auto r = item.reader->Get(item.id, i);
            assert(r);
            if (r)
                item.seq = i.seq;
            }
        }
    

    const std::unordered_set<Seq::Id>& IdsInFile(const std::string &fname) const;

protected:
    std::string DetectFileType(const std::string &fname);
    Seq::Id Insert(std::string &&name, std::string &&seq);
    Seq::Id Insert(const SeqReader::Item &item, SeqReader *reader, int mode);

protected:
    std::mutex mutex_;
    std::vector<std::string> names_;
    std::unordered_map<std::string, Seq::Id> names_to_ids_;

    mutable std::vector<Item> items_;
    std::unordered_map<std::string, std::unordered_set<Seq::Id>> ids_in_file_;
    std::vector<SeqReader*> readers_;
};

#endif // FSA_READ_STORE_HPP
