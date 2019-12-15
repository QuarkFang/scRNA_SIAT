from sqlalchemy import create_engine, Column, String, Integer
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from csv import reader
import os


Base = declarative_base()


def to_dict(self):
    return {c.name: getattr(self, c.name, None) for c in self.__table__.columns}


Base.to_dict = to_dict


class Marker(Base):
    __tablename__ = 'cell_marker'

    record_ID = Column(Integer, primary_key=True)
    species_type = Column(String(128))
    tissue_type = Column(String(128))
    uberon_ontology_ID = Column(String(128))
    cancer_type = Column(String(128))
    cell_type = Column(String(128))
    cell_name = Column(String(128))
    cell_ontology_ID = Column(String(128))
    cell_marker = Column(String(256))
    PMID = Column(String(128))


engine = create_engine('sqlite:///db/marker.db')
Base.metadata.create_all(engine)
DBSession = sessionmaker(bind=engine)


def add(record_ID: int, species_type: str, tissue_type: str, uberon_ontology_ID: str,
        cancer_type: str, cell_type: str, cell_name: str, cell_ontology_ID: str,
        cell_marker: str, PMID: str):

    session = DBSession()

    new = Marker(record_ID=record_ID, species_type=species_type, tissue_type=tissue_type,
                 uberon_ontology_ID=uberon_ontology_ID, cancer_type=cancer_type,
                 cell_type=cell_type, cell_name=cell_name, cell_ontology_ID=cell_ontology_ID,
                 cell_marker=cell_marker, PMID=PMID)
    session.add(new)
    session.commit()
    session.close()


def delete(record_ID: int):
    session = DBSession()
    marker = session.query(Marker).get(record_ID)
    session.delete(marker)
    session.commit()
    session.close()


def query(cell_marker):
    session = DBSession()
    result = session.query(Marker).filter(Marker.cell_marker.like('%'+cell_marker+'%')).all()
    cells = [cell.to_dict() for cell in result]
    session.commit()
    session.close()
    return cells


def update_db(filename: str):
    if not os.path.exists(filename):
        raise ValueError('the file is not existed!')

    session = DBSession()

    file = open(filename, 'r')
    data = reader(file)
    for i, row in enumerate(data):
        if i is not 0:
            record_ID = i - 1
            species_type = row[0]
            tissue_type = row[1]
            uberon_ontology_ID = row[2]
            cancer_type = row[3]
            cell_type = row[4]
            cell_name = row[5]
            cell_ontology_ID = row[6]
            cell_marker = row[7]
            PMID = row[13]

            new = Marker(record_ID=record_ID, species_type=species_type, tissue_type=tissue_type,
                         uberon_ontology_ID=uberon_ontology_ID, cancer_type=cancer_type,
                         cell_type=cell_type, cell_name=cell_name, cell_ontology_ID=cell_ontology_ID,
                         cell_marker=cell_marker, PMID=PMID)
            session.add(new)

    session.commit()
    session.close()
