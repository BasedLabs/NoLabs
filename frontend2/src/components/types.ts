export const PdbViews = {
    default: { key: 'base', title: 'Default' },
    cartoon: { key: 'cartoon', title: 'Cartoon' },
    backbone: { key: 'backbone', title: 'Backbone' },
    ballsAndSticks: { key: 'ball+stick', title: 'Balls and sticks' },
    contact: { key: 'contact', title: 'Contact' },
    helixorient: { key: 'helixorient', title: 'Helixorient' },
    hyperball: { key: 'hyperball', title: 'Hyperball' },
    licorice: { key: 'licorice', title: 'Licorice' },
    ribbon: { key: 'ribbon', title: 'Ribbon' },
    rope: { key: 'rope', title: 'Rope' },
    surface: { key: 'surface', title: 'Surface' },
    spacefill: { key: 'spacefill', title: 'Spacefill' },
    unitcell: { key: 'unitcell', title: 'Unitcell' }
};

export type ExperimentListItem = {
    id: string;
    name: string;
} | null;

export type Timeline = {
    message: string | null | undefined;
    error: string | null | undefined;
    createdAt: string | undefined;
}
