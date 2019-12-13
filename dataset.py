from sqlalchemy import create_engine, Column, String, Integer
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

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
    cell_marker = Column(String(128))
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
    result = session.query(Marker).filter(Marker.cell_marker == cell_marker).all()
    cells = [cell.to_dict() for cell in result]
    session.commit()
    session.close()
    return cells
