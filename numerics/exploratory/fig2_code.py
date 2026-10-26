import json
rows=json.load(open('fig2_rows.json'))
def fmt(v):
    s=('%.4f'%v).rstrip('0').rstrip('.')
    return s if s else '0'
rs=','.join("{phi:%s, slug:'%s', share:%s}"%(fmt(r['phi']),r['slug'],fmt(r['share'])) for r in rows)
designs="[{slug:'eq', design:'Equal loadings .70', short:'Equal loadings'},{slug:'mir', design:'Mirrored loadings .60 to .80', short:'Mirrored'},{slug:'mild', design:'Unequal blocks, .68 vs .72', short:'Mild'},{slug:'strong', design:'Unequal blocks, .60 vs .80', short:'Strong'}]"
marks="[{slug:'eq', phi:0.4089},{slug:'mir', phi:0.4411}]"
loc=json.load(open('fig2_local.json'))
ls=','.join("{phi:%s, slug:'mild', share:%s}"%(fmt(r['phi']),fmt(r['share'])) for r in loc)
jsx = (
"export default () => <claude.Visualize data-claude-component='tilt-beyond-exchangeability' sources={{rows:{kind:'data', data:[" + rs + "]}, designs:{kind:'data', data:" + designs + "}, marks:{kind:'data', data:" + marks + "}, local:{kind:'data', data:[" + ls + "]}}}>{({rows, designs, marks, local, datum}) => { "
"const left = 64, xr = 700, y0 = 100, yb = 410; "
"const x = f => left + f / 0.9 * (xr - left); const y = v => yb - (v - 0.5) / 0.5 * (yb - y0); "
"const dec = v => v.toFixed(3).replace('0.', '.'); const dec2 = v => v.toFixed(2).replace('0.', '.'); "
"const col = s => (s === 'mild' || s === 'strong') ? 'var(--cds-chart-categorical-1)' : 'var(--cds-chart-muted)'; "
"const op = s => s === 'mild' ? 0.5 : 1; const dash = s => s === 'mir' ? '6 4' : 'none'; const wid = s => s === 'strong' ? 2.5 : 2; "
"const name = s => designs.find(d => d.slug === s); const fold = local.at(-1); "
"return <svg viewBox='0 0 760 476' role='img' aria-label='Unequal blocks lean toward the stronger trait at every correlation' fontSize='12'>"
"<text data-claude-text-id='title' x='20' y='30' fontSize='16' fontWeight='600' fill='var(--cds-text-primary)'>Figure 2. Unequal blocks lean toward the stronger trait at every correlation</text>"
"<text data-claude-text-id='subtitle' x='20' y='54' fill='var(--cds-text-secondary)'>Share of the followed trait in the one-factor ML composite, 4 items per trait. Dotted: a local, non-global minimum.</text>"
"<text data-claude-text-id='subtitle-2' x='20' y='72' fill='var(--cds-text-secondary)'>.5 = equal composite, 1 = one trait only. Weights solve l_A = λ_A(w_A + φw_B) and l_B = λ_B(φw_A + w_B).</text>"
"<g data-claude-anchor='y-axis'><line x1={left} x2={left} y1={y0} y2={yb} stroke='var(--cds-chart-axis)'/>{[0.5, 0.6, 0.7, 0.8, 0.9, 1].map(t => <g key={t}><line x1={left - 4} x2={left} y1={y(t)} y2={y(t)} stroke='var(--cds-chart-axis)'/><text x={left - 8} y={y(t) + 4} textAnchor='end' fill='var(--cds-text-secondary)'>{t === 1 ? '1' : t.toFixed(1).replace('0.', '.')}</text></g>)}<text data-claude-text-id='y-title' x='16' y='255' transform='rotate(-90 16 255)' textAnchor='middle' fill='var(--cds-text-secondary)'>Share of the followed trait</text></g>"
"<g data-claude-anchor='x-axis'><line x1={left} x2={xr} y1={yb} y2={yb} stroke='var(--cds-chart-axis)'/>{[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].map(t => <text key={t} x={x(t)} y={yb + 18} textAnchor='middle' fill='var(--cds-text-secondary)'>{t === 0 ? '0' : t.toFixed(1).replace('0.', '.')}</text>)}<text data-claude-text-id='x-title' x='382' y='462' textAnchor='middle' fill='var(--cds-text-secondary)'>Trait correlation φ</text></g>"
"<g data-claude-anchor='thresholds'>{marks.map(mk => <g key={mk.slug}>{mk.slug === 'eq' ? <line x1={x(mk.phi)} x2={x(mk.phi)} y1={y0} y2={yb} stroke='var(--cds-chart-reference)' strokeDasharray='4 4' {...datum(mk, 'phi')}/> : null}<circle cx={x(mk.phi)} cy={y(0.5)} r='4' fill='var(--cds-chart-reference)' {...datum(mk, 'phi')}><title>{`${name(mk.slug).design}: threshold φ = ${dec(mk.phi)}`}</title></circle><text x={x(mk.phi) + 6} y={mk.slug === 'eq' ? y0 + 12 : y(0.5) - 9} fill='var(--cds-text-secondary)' {...datum(mk, 'phi')}>{mk.slug === 'eq' ? `${name(mk.slug).short}: φ* = ${dec(mk.phi)} (Theorem 1)` : `${name(mk.slug).short}: ${dec(mk.phi)}`}</text></g>)}</g>"
"{designs.map(d => { const pts = rows.filter(r => r.slug === d.slug); return <g key={d.slug} data-claude-anchor={`series-${d.slug}`}><path d={pts.map((r, i) => `${i === 0 ? 'M' : 'L'}${x(r.phi)},${y(r.share)}`).join(' ')} fill='none' stroke={col(d.slug)} strokeOpacity={op(d.slug)} strokeWidth={wid(d.slug)} strokeDasharray={dash(d.slug)}/>{pts.map(r => <circle key={r.phi} cx={x(r.phi)} cy={y(r.share)} r='1.8' fill={col(d.slug)} fillOpacity={op(d.slug)} {...datum(r, 'share')}><title>{`${d.design}, φ = ${dec(r.phi)}: share ${dec(r.share)}`}</title></circle>)}</g>; })}"
"<g data-claude-anchor='series-mild-local'><path d={local.map((r, i) => `${i === 0 ? 'M' : 'L'}${x(r.phi)},${y(r.share)}`).join(' ')} fill='none' stroke='var(--cds-chart-categorical-1)' strokeOpacity='0.5' strokeWidth='2' strokeDasharray='2 3'/>{local.map(r => <circle key={r.phi} cx={x(r.phi)} cy={y(r.share)} r='1.5' fill='var(--cds-chart-categorical-1)' fillOpacity='0.5' {...datum(r, 'share')}><title>{`${name(r.slug).design}, local minimum on the weaker trait, φ = ${dec(r.phi)}: share ${dec(r.share)}`}</title></circle>)}<circle cx={x(fold.phi)} cy={y(fold.share)} r='4' fill='none' stroke='var(--cds-chart-categorical-1)' strokeWidth='1.5' {...datum(fold, 'phi')}><title>{`Fold: the local minimum disappears at φ = ${dec(fold.phi)}`}</title></circle><text data-claude-text-id='local-label' x={x(0.06)} y={y(fold.share) + 22} fill='var(--cds-text-secondary)'>Local minimum on the weaker trait</text><text x={x(0.06)} y={y(fold.share) + 38} fill='var(--cds-text-secondary)' {...datum(fold, 'phi')}>{`${name(fold.slug).short} design, fold at φ = ${dec2(fold.phi)}`}</text></g>"
"<g data-claude-anchor='end-values'>{rows.filter(r => (r.slug === 'mild' || r.slug === 'strong') && Math.abs(r.phi - 0.9) < 1e-9).map(r => <text key={r.slug} x={x(r.phi) + 6} y={y(r.share) + 4} fill='var(--cds-text-primary)' {...datum(r, 'share')}>{dec2(r.share)}</text>)}</g>"
"<g data-claude-anchor='legend'>{designs.map((d, i) => <g key={d.slug}><line x1={left + 18} x2={left + 44} y1={yb - 76 + i * 18} y2={yb - 76 + i * 18} stroke={col(d.slug)} strokeOpacity={op(d.slug)} strokeWidth={wid(d.slug)} strokeDasharray={dash(d.slug)}/><text x={left + 52} y={yb - 72 + i * 18} fill='var(--cds-text-primary)'>{d.design}</text></g>)}</g>"
"</svg>; }}</claude.Visualize>;"
)
open('fig2_code.jsx','w').write(jsx)
print(len(jsx)); print(jsx[-1500:])
